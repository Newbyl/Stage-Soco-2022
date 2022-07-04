#include "Cplate.h"

// TODO : Compile with /O2 /fp:fast compilation option for windows (MSCV compiler).
// TODO : Compile with -O3 -ffast-math flags for linux / mac (g++ / clang++ compiler).
// These options / flags are used for performance gain (known bug : if the values are completely false, remove the flags, compile, put the flags again and then recompile.
// This bug was experienced on macOS Monterey with clang++ compiler and i dont know why this happened).

Cplate::Cplate()
{
	// Initialize Cplate's members

	opened = false;
	calculationDone = false;

	resolution = 0.1;

	numberOfElements = 0;

	centerApertureX = 0;
	centerApertureY = 0;
	centerApertureZ = 0;

	elements.coordinates.x = NULL;
	elements.coordinates.y = NULL;
	elements.coordinates.z = NULL;

	material.velocity = 3900.0;

	coupling.velocity = 1500.0;
	coupling.height = 0.0;

	positionType = FLOW_POSITION_TYPE::PATH;

	numberOfTargets = 0;

	targets.tilts = NULL;
	targets.skews = NULL;
	targets.positions = NULL;
	laws = NULL;
	paths = NULL;
}

Cplate::~Cplate()
{
	Close();
}

int Cplate::Open()
{
	if (opened)
		return PLUGIN_ALREADY_OPEN;

	opened = true;
	return PLUGIN_NO_ERROR;
}

int Cplate::Close()
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	opened = false;

	// Release allocated memories...

	if (targets.tilts != NULL)
		free(targets.tilts);
	if (targets.skews != NULL)
		free(targets.skews);
	if (targets.positions != NULL)
		free(targets.positions);

	for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
	{
		if (laws != NULL)
			if (laws[iTarget].delays != NULL)
				free(laws[iTarget].delays);
	}
	if (laws != NULL)
		free(laws);
	if (paths != NULL)
		free(paths);

	if (elements.coordinates.x != NULL)
		free(elements.coordinates.x);
	if (elements.coordinates.y != NULL)
		free(elements.coordinates.y);
	if (elements.coordinates.z != NULL)
		free(elements.coordinates.z);

	return PLUGIN_NO_ERROR;
}

int Cplate::Set(const char* param_name, int unit, int* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "Flow.PositionType") == 0)
	{
		positionType = (FLOW_POSITION_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}

int Cplate::Set(const char* param_name, int unit, double* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "Material.Velocity") == 0)
	{
		material.velocity = Unit::ChangeUnit(*value, unit, UNIT_mps);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Coupling.Velocity") == 0)
	{
		coupling.velocity = Unit::ChangeUnit(*value, unit, UNIT_mps);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Coupling.Height") == 0)
	{
		coupling.height = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "center.aperture.x") == 0)
	{
		centerApertureX = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "center.aperture.y") == 0)
	{
		centerApertureY = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "center.aperture.z") == 0)
	{
		centerApertureZ = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cplate::Set(const char* param_name, int unit, char* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Set(const char* param_name, int unit, int* dim1, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Set(const char* param_name, int unit, int* dim1, double* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strnicmp(param_name, "Elements.Coordinates", strlen("Elements.Coordinates")) == 0)
	{
		if (numberOfElements > 0 && *dim1 != numberOfElements) // Do we have to resize arrays of coordinates?
		{													   // Yes. Release previous arrays.
			free(elements.coordinates.x);
			free(elements.coordinates.y);
			free(elements.coordinates.z);
			elements.coordinates.x = NULL;
			elements.coordinates.y = NULL;
			elements.coordinates.z = NULL;
			numberOfElements = *dim1;
		}
		if (numberOfElements == 0)
			numberOfElements = *dim1;
	}

	if (strcmpi(param_name, "Elements.Coordinates.x") == 0)
	{
		if (elements.coordinates.x == NULL)
			elements.coordinates.x = (double*)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.x[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Elements.Coordinates.y") == 0)
	{
		if (elements.coordinates.y == NULL)
			elements.coordinates.y = (double*)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.y[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Elements.Coordinates.z") == 0)
	{
		if (elements.coordinates.z == NULL)
			elements.coordinates.z = (double*)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.z[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strnicmp(param_name, "Targets", strlen("Targets")) == 0)
	{
		if (numberOfTargets > 0 && *dim1 != numberOfTargets) // Do we have to resize arrays of targets?
		{													 // Yes. Release previous arrays.
			free(targets.tilts);
			free(targets.skews);
			free(targets.positions);
			free(laws);
			free(paths);
			targets.tilts = NULL;
			targets.skews = NULL;
			targets.positions = NULL;
			laws = NULL;
			paths = NULL;
			numberOfTargets = *dim1;
		}
		if (numberOfTargets == 0)
			numberOfTargets = *dim1;
	}

	if (strcmpi(param_name, "Targets.Tilts") == 0)
	{
		if (targets.tilts == NULL)
			targets.tilts = (double*)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.tilts[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_rad);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Targets.Skews") == 0)
	{
		if (targets.skews == NULL)
			targets.skews = (double*)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.skews[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_rad);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Targets.Positions") == 0)
	{
		if (targets.positions == NULL)
			targets.positions = (double*)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.positions[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	calculationDone = false;
	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cplate::Set(const char* param_name, int unit, int* dim1, int* dim2, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Set(const char* param_name, int unit, int* dim1, int* dim2, double* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Set(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Set(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, double* value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cplate::Get(const char* param_name, int unit, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Get(const char* param_name, int unit, double* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	// TODO Ceci est un exemple pour v�rifier la fonction Get! Nous n'avons appriori pas besoin de relire les valeurs des param�tres!
	if (strcmpi(param_name, "Material.Velocity") == 0)
	{
		*value = Unit::ChangeUnit(material.velocity, UNIT_mps, unit);
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cplate::Get(const char* param_name, int unit, int* dim1, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Get(const char* param_name, int unit, int* dim1, double* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Get(const char* param_name, int unit, int* dim1, char* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "version") == 0)
	{
		if (dim1 == NULL)
			return PLUGIN_INVALID_ARGUMENT;
		if (value == NULL || *dim1 < strlen("1.0.0.1"))
			*dim1 = strlen("1.0.0.1");
		else
			strcpy(value, "1.0.0.1");
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cplate::Get(const char* param_name, int unit, int* dim1, int* dim2, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Get(const char* param_name, int unit, int* dim1, int* dim2, double* value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "Laws") == 0)
	{
		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < numberOfElements)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iElement = 0; iElement < numberOfElements; iElement++)
			{
				*(value + iTarget * numberOfElements + iElement) = Unit::ChangeUnit(laws[iTarget].delays[iElement], UNIT_us, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Paths.x") == 0)
	{
		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < 3)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < 3; iPoint++)
			{
				*(value + iTarget * 3 + iPoint) = Unit::ChangeUnit(paths[iTarget].x[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}
	if (strcmpi(param_name, "Paths.y") == 0)
	{
		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < 3)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < 3; iPoint++)
			{
				*(value + iTarget * 3 + iPoint) = Unit::ChangeUnit(paths[iTarget].y[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}
	if (strcmpi(param_name, "Paths.z") == 0)
	{
		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < 3)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < 3; iPoint++)
			{
				*(value + iTarget * 3 + iPoint) = Unit::ChangeUnit(paths[iTarget].z[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cplate::Get(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, int* value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cplate::Get(const char* param_name, int unit, int* dim1, int* dim2, int* dim3, double* value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cplate::ExecSync(const char* action)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(action, "Calculate") == 0)
	{
		// Check if the calculation is possible
		if (calculationDone)
			return PLUGIN_NO_ERROR; // Calculation is already done!
		if (numberOfElements == 0 || elements.coordinates.x == NULL || elements.coordinates.y == NULL || elements.coordinates.z == NULL)
			return PLUGIN_PARMETER_UNDEFINED;
		if (numberOfTargets == 0 || targets.tilts == NULL || targets.skews == NULL || targets.positions == NULL)
			return PLUGIN_PARMETER_UNDEFINED;

		// Allocate resources for laws and paths
		laws = (Cplate::PLAW)malloc(numberOfTargets * sizeof(Cplate::LAW));
		paths = (Cplate::PPATH)malloc(numberOfTargets * sizeof(Cplate::PATH));
		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			laws[iTarget].delays = (double*)malloc(numberOfElements * sizeof(double));
		}

		// Calculates laws and paths
		Calculate();

		calculationDone = true;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}

int Cplate::ExecAsync(const char* action) { return PLUGIN_UNKNOWN_PARAMETER; }


// Get the maximum element of an array.
double Cplate::maxArray(double* array, int size)
{
	double max = array[0];

	for (int i = 0; i < size; i++)
	{
		if (max <= array[i])
		{
			max = array[i];
		}
	}

	return max;
}


// Get the minimum element of an array.
double Cplate::minArray(double* array, int size)
{
	double min = array[0];

	for (int i = 0; i < size; i++)
	{
		if (min > array[i])
		{
			min = array[i];
		}
	}

	return min;
}


// Get the index of the minimum element of an array.
int Cplate::minArrayIndex(const double* array, int size)
{
	double min = array[0];
	int index = 0;

	for (int i = 0; i < size; i++)
	{
		if (min > array[i])
		{
			min = array[i];
			index = i;
		}
	}

	return index;
}


int Cplate::Calculate()
{
	// This calculation is based on a probe centered on X=Y=Z=0.0
	if (coupling.height > 0.0)
	{
		for (int i = 0; i < numberOfElements; i++)
		{
			elements.coordinates.y[i] = coupling.height - elements.coordinates.y[i];
		}
	}

	// Here we get the max and min of x,y,z coordinates.
	double maxYProbe = maxArray(elements.coordinates.y, numberOfElements);
	double minYProbe = minArray(elements.coordinates.y, numberOfElements);
	double maxXProbe = maxArray(elements.coordinates.x, numberOfElements);
	double minXProbe = minArray(elements.coordinates.x, numberOfElements);
	double maxZProbe = maxArray(elements.coordinates.z, numberOfElements);
	double minZProbe = minArray(elements.coordinates.z, numberOfElements);

	// Memory allocation for positionProjection array (temporary array to not modify the values int targets.position array).
	double* positionProjection = (double*)malloc(numberOfTargets * sizeof(double));
	// Memory allocation for positionProjection2 array for all the focalisation point for remarkable points.
	double* positionProjection2 = (double*)malloc(numberOfTargets * sizeof(double));

	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
		if (positionType == FLOW_POSITION_TYPE::DEPTH)
		{
			positionProjection2[iLaw] = targets.positions[iLaw];
			positionProjection[iLaw] = targets.positions[iLaw];
		}
		else {
			positionProjection2[iLaw] = targets.positions[iLaw] * cos(targets.tilts[iLaw]);
			positionProjection[iLaw] = targets.positions[iLaw] / cos(targets.tilts[iLaw]);
		}
	}

	// This condition handle the case when the distance between the probe and the plate equals to 0.
	if (coupling.height == 0)
	{
		// Array of all distances between probe element and focal point.
		double* distances;

		// For loop that iterate the number of laws that we want.
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			// Memory allocation for distances array.
			distances = (double*)malloc(numberOfElements * sizeof(double));

			// For loop that compute the distance between probe element and focal point and push it into distances array.
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				distances[iElem] = sqrt(pow((cos(targets.tilts[iLaw]) * positionProjection[iLaw]) - elements.coordinates.y[iElem], 2.0) + pow((cos(targets.skews[iLaw]) * (sin(targets.tilts[iLaw]) * positionProjection[iLaw])) - elements.coordinates.x[iElem] + (((maxXProbe - minXProbe) / 2) + minXProbe), 2.0) + pow(((sin(targets.tilts[iLaw]) * positionProjection[iLaw]) * sin(targets.skews[iLaw])) - elements.coordinates.z[iElem], 2.0));
			}

			// Maximum value of distances array.
			double maxDistances = maxArray(distances, numberOfElements);

			// Compute the delay for each element of the probe, convert it to micro seconds and push it into laws array.
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				laws[iLaw].delays[iElem] = ((maxDistances - distances[iElem]) / material.velocity) * 1000;
			}

			// Computation of the remarkable points.
			paths[iLaw].x[0] = centerApertureX;
			paths[iLaw].x[1] = centerApertureX;
			paths[iLaw].x[2] = (sin(targets.tilts[iLaw]) * targets.positions[iLaw]) * cos(targets.skews[iLaw]) + paths[iLaw].x[0];

			paths[iLaw].y[0] = -centerApertureY;
			paths[iLaw].y[1] = -centerApertureY;
			paths[iLaw].y[2] = positionProjection2[iLaw] + paths[iLaw].y[0];

			paths[iLaw].z[0] = centerApertureZ;
			paths[iLaw].z[1] = centerApertureZ;
			paths[iLaw].z[2] = paths[iLaw].x[2] * (tan(targets.skews[iLaw]));
			

			// Release of the memory taken by distances array.
			free(distances);
		}
	}

	// This else handle the case when coupling height do not equals to 0.
	else
	{
		// Memory allocation for x,y,z coordinates of flat bottom holes.
		double* xFbhP = (double*)malloc(numberOfTargets * sizeof(double));
		double* yFbhP = (double*)malloc(numberOfTargets * sizeof(double));
		double* zFbhP = (double*)malloc(numberOfTargets * sizeof(double));

		// For loop that compute x,y,z coordinates of flat bottom holes and push them into their respectives array.
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			double yFhb = (cos(targets.tilts[iLaw]) * targets.positions[iLaw]) + coupling.height;

			double xFhb = ((tan(asin((sin(targets.tilts[iLaw]) * coupling.velocity) / material.velocity)) * coupling.height) + (sin(targets.tilts[iLaw]) * targets.positions[iLaw])) * cos(targets.skews[iLaw])
				+ (((maxXProbe - minXProbe) / 2) + minXProbe);

			double zFhb = ((tan(asin((sin(targets.tilts[iLaw]) * coupling.velocity) / material.velocity)) * coupling.height) + (sin(targets.tilts[iLaw]) * targets.positions[iLaw])) * sin(targets.skews[iLaw]);

			xFbhP[iLaw] = xFhb;
			yFbhP[iLaw] = yFhb;
			zFbhP[iLaw] = zFhb;
		}

		// Memory allocation array1, array2, array3, array4 are where we are going to put all min / max boundaries for interest zones on the interface.
		double* array1 = (double*)malloc(numberOfTargets * sizeof(double));
		double* array2 = (double*)malloc(numberOfTargets * sizeof(double));
		double* array3 = (double*)malloc(numberOfTargets * sizeof(double));
		double* array4 = (double*)malloc(numberOfTargets * sizeof(double));

		// Creation of counters to know the length of aboves array.
		int cptAr1 = 0;
		int cptAr2 = 0;
		int cptAr3 = 0;
		int cptAr4 = 0;

		// For loop that determines min / max boundaries for the interest zone on the interface.
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			if (maxXProbe < xFbhP[iLaw])
			{
				array1[iLaw] = ((coupling.height - minYProbe) * ((xFbhP[iLaw] - maxXProbe) / (yFbhP[iLaw] - minYProbe))) + maxXProbe;
				cptAr1++;
			}
			else
			{
				array1[iLaw] = maxXProbe;
				cptAr1++;
			}

			if (xFbhP[iLaw] < minXProbe)
			{
				array2[iLaw] = (minXProbe + xFbhP[iLaw]) / (yFbhP[iLaw] - minYProbe) * (coupling.height - minYProbe) + minXProbe;
				cptAr2++;
			}
			else
			{
				array2[iLaw] = minXProbe;
				cptAr2++;
			}

			if (zFbhP[iLaw] > maxZProbe)
			{
				array3[iLaw] = (maxZProbe + ((coupling.height - minYProbe) * ((zFbhP[iLaw] - maxZProbe) / (yFbhP[iLaw] - minYProbe))));
				cptAr3++;
			}
			else
			{
				array3[iLaw] = maxZProbe;
				cptAr3++;
			}

			if (zFbhP[iLaw] < minZProbe)
			{
				array4[iLaw] = (((zFbhP[iLaw] + minZProbe) / (yFbhP[iLaw] - minYProbe)) * (coupling.height - minYProbe)) + minZProbe;
				cptAr4++;
			}
			else
			{
				array4[iLaw] = minZProbe;
				cptAr4++;
			}
		}

		double interestZoneSize1 = maxArray(array1, cptAr1) - minArray(array2, cptAr2); // Size of the interest zone on the XY axis.
		double interestZoneSize2 = maxArray(array3, cptAr3) - minArray(array4, cptAr4); // Size of the interest zone on the YZ axis.


		// Release of the memory taken by array1 and array3.
		free(array1);
		free(array3);

		std::vector<double> arrayCouplingHeight; // Array where we are going to put the coupling height before duplicating it into xIntP.
		std::vector<double> arrayPreYInterfaceP; // Array that contain Y coodinates of interest points before duplicating it into yIntP.

		// For loop that compute Y interface points coordinates and X interface points coordinates and push them into arrays (Y interface points here are contained in arrayPreInterfaceP
		// and X interface points in arrayCouplingHeight, explications in the comment below.
		for (int i = 0; i < (int)((interestZoneSize1 / resolution) + 1); i++)
		{
			// Initially the algorithm was calculating interface points on YX axis but by convention we take XY that's why the names of arrays are inverted here.
			arrayCouplingHeight.push_back((resolution * i) + minArray(array2, cptAr2));
			arrayPreYInterfaceP.push_back(coupling.height);
		}


		// Release of the memory taken by array2.
		free(array2);

		// Memory allocation for x,y,z coordinates of interface points array.
		std::vector<double> xIntP;
		std::vector<double> yIntP;
		std::vector<double> zIntP;

		// For loop that duplicate X,Y coordinates and compute Z interface point coordinate for matrix probe.
		for (int i = 0; i < (int)((interestZoneSize2 / resolution) + 1); i++)
		{
			for (int j = 0; j < arrayCouplingHeight.size(); j++)
			{
				xIntP.push_back(arrayCouplingHeight[j]);
				yIntP.push_back(arrayPreYInterfaceP[j]);
			}

			for (int j = 0; j < arrayPreYInterfaceP.size(); j++)
			{
				zIntP.push_back((i * resolution) + minArray(array4, cptAr4));
			}
		}

		// Release of the memory taken by array4.
		free(array4);

		// This if not define is the version of the algorithm non-optimized with dichotomy.
		// If you want to use this one, put in comment "#define _OPTIMIZATION" in "Cplate.h".
		#ifndef _OPTIMIZATION
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			double* timeFbhInt = (double*)malloc((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the interest point to the flat bottom hole

			// For loop that compute the time for the US to go from the interest point to the flat bottom hole
			for (int j = 0; j < (int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1); j++)
			{
				timeFbhInt[j] = (sqrt(pow(xFbhP[iLaw] - xIntP[j], 2.0) + pow(yFbhP[iLaw] - yIntP[j], 2.0) + pow(zFbhP[iLaw] - zIntP[j], 2.0))) / (material.velocity / 1000);
			}

			double* delayLaw = (double*)malloc(numberOfElements * sizeof(double));																  // Array that contain delay law for each probe
			double* addTimeVector = (double*)malloc((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the probe element to the flat bottom hole

			int minIndex;

			// For loop that compute the time for the US to go from the probe element to the interest point
			for (int probeElem = 0; probeElem < numberOfElements; probeElem++)
			{
				double* timeProbeInt = (double*)malloc((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the probe element to the interest point

				// For loop that compute the time for the US to go from the probe element to the interest point
				for (int j = 0; j < (int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1); j++)
				{
					timeProbeInt[j] = (sqrt(pow(elements.coordinates.x[probeElem] - xIntP[j], 2.0) + pow(elements.coordinates.y[probeElem] - yIntP[j], 2.0) + pow(elements.coordinates.z[probeElem] - zIntP[j], 2.0))) / (coupling.velocity / 1000);
				}

				// For loop that add timeFbhInt and timeProbeInt
				for (int l = 0; l < (int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1); l++)
				{
					addTimeVector[l] = timeFbhInt[l] + timeProbeInt[l];
				}

				free(timeProbeInt);

				delayLaw[probeElem] = minArray(addTimeVector, (int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1));
			}

			minIndex = minArrayIndex(addTimeVector, (int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1));

			free(timeFbhInt);
			free(addTimeVector);

			// Release the memory taken by timeFbhInt.

			// Maximum element of delayLaw array.
			double maxDelayLaw = maxArray(delayLaw, numberOfElements);
			// For loop that push the delay law for each element of the probe
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{

				laws[iLaw].delays[iElem] = maxDelayLaw - delayLaw[iElem];
			}

			// Release of the memory taken by delayLaw array.
			free(delayLaw);

			double focusPointDistance = (sqrt(pow(xFbhP[iLaw] - xIntP[minIndex], 2.0) + pow(yFbhP[iLaw] - yIntP[minIndex], 2.0) + pow(zFbhP[iLaw] - zIntP[minIndex], 2.0)));

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = centerApertureX;
				paths[iLaw].x[1] = tan(asin((coupling.velocity * sin(targets.tilts[iLaw])) / material.velocity)) * coupling.height + paths[iLaw].x[0];
				paths[iLaw].x[2] = (sin(targets.tilts[iLaw]) * targets.positions[iLaw]) * cos(targets.skews[iLaw]) + paths[iLaw].x[0];

				paths[iLaw].y[0] = -centerApertureY;
				paths[iLaw].y[1] = coupling.height + paths[iLaw].y[0];
				paths[iLaw].y[2] = positionProjection2[iLaw];

				paths[iLaw].z[0] = centerApertureZ;
				paths[iLaw].z[1] = sin(targets.skews[iLaw]) * paths[iLaw].x[1];
				paths[iLaw].z[2] = paths[iLaw].x[2] * (tan(targets.skews[iLaw]));
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = centerApertureX;
				paths[iLaw].x[1] = tan(asin((coupling.velocity * sin(targets.tilts[iLaw])) / material.velocity)) * coupling.height + paths[iLaw].x[0];
				paths[iLaw].x[2] = (sin(targets.tilts[iLaw]) * targets.positions[iLaw]) * cos(targets.skews[iLaw]) + paths[iLaw].x[0];

				paths[iLaw].y[0] = -centerApertureY;
				paths[iLaw].y[1] = coupling.height + paths[iLaw].y[0];
				paths[iLaw].y[2] = positionProjection2[iLaw];

				paths[iLaw].z[0] = centerApertureZ;
				paths[iLaw].z[1] = sin(targets.skews[iLaw]) * paths[iLaw].x[1];
				paths[iLaw].z[2] = paths[iLaw].x[2] * (tan(targets.skews[iLaw]));
			}
		}
		#endif

		// This if define is the version of the algorithm optimized with dichotomy.
		// If you want to use this one, delete the double slash before "define _OPTIMIZATION".
		#ifdef _OPTIMIZATION
		// For loop that iterates number of law times.
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			// Array where we gonna put all our delay for each probe.
			double* delay = (double*)malloc(numberOfElements * sizeof(double));

			// Array that contain the time for the US to go from the interface point to the flat bottom hole.
			double* timeFbhInt = (double*)malloc((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1) * sizeof(double));

			// For loop that compute the time for the US to go from the interface point to the flat bottom hole.
			for (int j = 0; j < ((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1)); j++)
			{
				timeFbhInt[j] = (sqrt(pow(xFbhP[iLaw] - xIntP[j], 2.0) + pow(yFbhP[iLaw] - yIntP[j], 2.0) + pow(zFbhP[iLaw] - zIntP[j], 2.0))) / (material.velocity / 1000);
			}

			int nbGroupInt = (int)((interestZoneSize2 / resolution) + 1);

			double addTimeElemIntFbh = 0; // Addition of the time taken by the US between the probe element and the interface and between the interface to the flat bottom hole.

			// For loop that initialize each value of delay array to infinity.
			for (int i = 0; i < numberOfElements; i++)
			{
				delay[i] = INFINITY;
			}

			// This variable gives us the offset between each parabola of values of interface points.
			double offset = xIntP.size() / nbGroupInt;

			// For loop that iterate number of "groups" of interface times to get the minimum value between all interface points.
			for (int i = 0; i < nbGroupInt; i++)
			{
				// For loop that compute the time for the US to go from the probe element to the interface point.
				for (int probeElem = 0; probeElem < numberOfElements; probeElem++)
				{
					double* timeProbeInt = (double*)malloc((int)((interestZoneSize2 / resolution) + 1) * (int)((interestZoneSize1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the probe element to each interface point.

					// Computation of the time between the probe element and the first interface point.
					timeProbeInt[0] = (sqrt(pow(elements.coordinates.x[probeElem] - xIntP[0], 2.0) + pow(elements.coordinates.y[probeElem] - yIntP[0], 2.0) + pow(elements.coordinates.z[probeElem] - zIntP[0], 2.0))) / (coupling.velocity / 1000);

					// Indexes and length for the while loop below.
					int start = (offset * i);
					int end = (offset * (i + 1)) - 1;
					int mid = (int)((end - start) / 2);
					int length = 0;

					// While loop the find the minimum distance by dichotomy.
					while (length != -1)
					{
						length = end - start;
						mid = start + (int)(length / 2);

						if (mid == 0 || mid == xIntP.size())
							break;

						// Check if the previous time taken by the US is lower than the current.
						if (timeFbhInt[mid - 1] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntP[mid - 1], 2.0) + pow(elements.coordinates.y[probeElem] - yIntP[mid - 1], 2.0) + pow(elements.coordinates.z[probeElem] - zIntP[mid - 1], 2.0))) / (coupling.velocity / 1000) <
							timeFbhInt[mid] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntP[mid], 2.0) + pow(elements.coordinates.y[probeElem] - yIntP[mid], 2.0) + pow(elements.coordinates.z[probeElem] - zIntP[mid], 2.0))) / (coupling.velocity / 1000))
						{
							end = mid - 1;
						}

						// If the previous time taken by the US is bigger than the current we enter in this case
						else
						{
							start = mid + 1;
						}
					}

					// The addition of the time taken by the US between the probe element and the interface and between the interface to the flat bottom hole
					// and convert it to mm per seconds.
					addTimeElemIntFbh = timeFbhInt[end] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntP[end], 2.0) + pow(elements.coordinates.y[probeElem] - yIntP[end], 2.0) + pow(elements.coordinates.z[probeElem] - zIntP[end], 2.0))) / (coupling.velocity / 1000);

					// Release of the memory taken by timeProbeInt.
					free(timeProbeInt);

					// Push addTimeElemIntFbh into delayLaw array.
					if (delay[probeElem] > addTimeElemIntFbh)
					{
						delay[probeElem] = addTimeElemIntFbh;
					}
				}
			}


			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = centerApertureX;
				paths[iLaw].x[1] = (tan(asin((coupling.velocity * sin(targets.tilts[iLaw])) / material.velocity)) * coupling.height + paths[iLaw].x[0]);
				paths[iLaw].x[2] = (sin(targets.tilts[iLaw]) * targets.positions[iLaw]) * cos(targets.skews[iLaw]) + paths[iLaw].x[0];

				paths[iLaw].y[0] = -centerApertureY;
				paths[iLaw].y[1] = coupling.height + paths[iLaw].y[0];
				paths[iLaw].y[2] = positionProjection2[iLaw];

				paths[iLaw].z[0] = centerApertureZ;
				paths[iLaw].z[1] = sin(targets.skews[iLaw]) * paths[iLaw].x[1];
				paths[iLaw].z[2] = paths[iLaw].x[2] * (tan(targets.skews[iLaw]));
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = centerApertureX;
				paths[iLaw].x[1] = (tan(asin((coupling.velocity * sin(targets.tilts[iLaw])) / material.velocity)) * coupling.height + paths[iLaw].x[0]);
				paths[iLaw].x[2] = ((sin(targets.tilts[iLaw]) * targets.positions[iLaw]) * cos(targets.skews[iLaw]) + paths[iLaw].x[0]);

				paths[iLaw].y[0] = -centerApertureY;
				paths[iLaw].y[1] = coupling.height + paths[iLaw].y[0];
				paths[iLaw].y[2] = positionProjection2[iLaw];

				paths[iLaw].z[0] = centerApertureZ;
				paths[iLaw].z[1] = sin(targets.skews[iLaw]) * paths[iLaw].x[1];
				paths[iLaw].z[2] = paths[iLaw].x[2] * (tan(targets.skews[iLaw]));
			}

			// Release the memory taken by timeFbhInt.
			free(timeFbhInt);

			// Maximum element of delay array.
			double maxDelay = maxArray(delay, numberOfElements);
			// For loop that push the delay law for each element of the probe
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				laws[iLaw].delays[iElem] = maxDelay - delay[iElem];
			}

			// Release of the memory taken by delayLaw array.
			free(delay);
		}
		#endif

		// Release of the memory taken by x/y/zFbhP
		free(xFbhP);
		free(yFbhP);
		free(zFbhP);
	}

	// Release of the memory taken by positionProjection array.
	free(positionProjection);
	free(positionProjection2);

	// Retreive original coorindates
	if (coupling.height > 0.0)
	{
		for (int i = 0; i < numberOfElements; i++)
		{
			elements.coordinates.y[i] = -1.0*elements.coordinates.y[i] + coupling.height;
		}
	}

	return PLUGIN_NO_ERROR;
}