#include "framework.h"
#include "Cbar.h"
#include "unit.h"
#include "error_codes.h"
#include <cstddef>
#include <string.h>
#include "framework.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#define _CRT_SECURE_NO_WARNINGS


#define M_PI 3.14159265358979323846264338327950288

#define _OPTIMIZATION

#define strcmpi(x, y) strcasecmp((x), (y))
#define _strnicmp(x, y, z) strncasecmp((x), (y), (z))

// Ajouter flag /fp:fast et /O2 a la compilation.

Cbar::Cbar()
{
	// Initialize Cplate's members

	opened = false;
	calculationDone = false;

	resolution = 0.1;
	barDiameter = 0;

	numberOfElements = 0;

	elements.coordinates.x = NULL;
	elements.coordinates.y = NULL;
	elements.coordinates.z = NULL;

	material.velocity = 3900.0;

	coupling.velocity = 1480.0;
	coupling.height = 0.0;

	angleType = ANGLE_TYPE::INCIDENT;
	defectType = DEFECT_TYPE::FBH;

	numberOfTargets = 0;

	targets.tilts = NULL;
	targets.positions = NULL;
	targets.notchesAngles = NULL;
	laws = NULL;
	paths = NULL;
}

Cbar::~Cbar()
{
	Close();
}

int Cbar::Open()
{
	if (opened)
		return PLUGIN_ALREADY_OPEN;

	opened = true;
	return PLUGIN_NO_ERROR;
}

int Cbar::Close()
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	opened = false;

	// Release allocated memories...

	if (targets.tilts != NULL)
		free(targets.tilts);
	if (targets.positions != NULL)
		free(targets.positions);
	if(targets.notchesAngles != NULL)
		free(targets.notchesAngles);

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

int Cbar::Set(const char *param_name, int unit, int *value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "Angle_Type") == 0)
	{
		angleType = (ANGLE_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Defect_Type") == 0)
	{
		defectType = (DEFECT_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}

int Cbar::Set(const char *param_name, int unit, double *value)
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

	if (strcmpi(param_name, "Angle_Type") == 0)
	{
		angleType = (ANGLE_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Defect_Type") == 0)
	{
		defectType = (DEFECT_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Resolution") == 0)
	{
		resolution = *value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "BarDiameter") == 0)
	{
		barDiameter = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cbar::Set(const char *param_name, int unit, char *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cbar::Set(const char *param_name, int unit, int *dim1, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cbar::Set(const char *param_name, int unit, int *dim1, double *value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (_strnicmp(param_name, "Elements.Coordinates", strlen("Elements.Coordinates")) == 0)
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
			elements.coordinates.x = (double *)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.x[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Elements.Coordinates.y") == 0)
	{
		if (elements.coordinates.y == NULL)
			elements.coordinates.y = (double *)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.y[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Elements.Coordinates.z") == 0)
	{
		if (elements.coordinates.z == NULL)
			elements.coordinates.z = (double *)malloc(*dim1 * sizeof(double));
		for (int iElement = 0; iElement < *dim1; iElement++)
			elements.coordinates.z[iElement] = Unit::ChangeUnit(value[iElement], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (_strnicmp(param_name, "Targets", strlen("Targets")) == 0)
	{
		if (numberOfTargets > 0 && *dim1 != numberOfTargets) // Do we have to resize arrays of targets?
		{													 // Yes. Release previous arrays.
			free(targets.tilts);
			free(targets.positions);
			free(targets.notchesAngles);
			free(laws);
			free(paths);
			targets.tilts = NULL;
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
			targets.tilts = (double *)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.tilts[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_deg);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Targets.NotchesAngles") == 0)
	{
		if (targets.notchesAngles == NULL)
			targets.notchesAngles = (double *)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.notchesAngles[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_deg);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Targets.Positions") == 0)
	{
		if (targets.positions == NULL)
			targets.positions = (double *)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.positions[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	calculationDone = false;
	return PLUGIN_UNKNOWN_PARAMETER;
}
int Cbar::Set(const char *param_name, int unit, int *dim1, int *dim2, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Set(const char *param_name, int unit, int *dim1, int *dim2, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Set(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Set(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cbar::Get(const char *param_name, int unit, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Get(const char *param_name, int unit, double *value)
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
int Cbar::Get(const char *param_name, int unit, int *dim1, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Get(const char *param_name, int unit, int *dim1, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Get(const char *param_name, int unit, int *dim1, char *value)
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
int Cbar::Get(const char *param_name, int unit, int *dim1, int *dim2, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Get(const char *param_name, int unit, int *dim1, int *dim2, double *value)
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
int Cbar::Get(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Cbar::Get(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Cbar::ExecSync(const char *action)
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
		if (numberOfTargets == 0 || targets.tilts == NULL || targets.positions == NULL)
			return PLUGIN_PARMETER_UNDEFINED;

		// Allocate resources for laws and paths
		laws = (Cbar::PLAW)malloc(numberOfTargets * sizeof(Cbar::LAW));
		paths = (Cbar::PPATH)malloc(numberOfTargets * sizeof(Cbar::PATH));
		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			laws[iTarget].delays = (double *)malloc(numberOfElements * sizeof(double));
		}

		// Calculates laws and paths
		Calculate();

		calculationDone = true;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}

int Cbar::ExecAsync(const char *action) { return PLUGIN_UNKNOWN_PARAMETER; }

/**
 * maxArray
 *
 * get the max element from a array
 *
 * @param array<double> array : The array you want to get the max element
 *
 * @return double max : The max element of the array
 **/
double Cbar::maxArray(double *array, int size)
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

/**
 * minArray
 *
 * get the min element from a array
 *
 * @param array<double> array : The array you want to get the min element
 *
 * @return double min : The min element of the array
 **/
double Cbar::minArray(double *array, int size)
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


std::vector<double*> Cbar::fbhBuilder(double barRadius)
{
	double* ai = (double*)malloc(numberOfTargets * sizeof(double));
	double* ar = (double*)malloc(numberOfTargets * sizeof(double));

	double* utAngle = (double*)malloc(numberOfTargets * sizeof(double));

	double* x = (double*)malloc(numberOfTargets * sizeof(double));
	double* y = (double*)malloc(numberOfTargets * sizeof(double));
	double* z = (double*)malloc(numberOfTargets * sizeof(double));

	double* xInt = (double*)malloc(numberOfTargets * sizeof(double));
	double* yInt = (double*)malloc(numberOfTargets * sizeof(double));
	double* zInt = (double*)malloc(numberOfTargets * sizeof(double));

	double minXprobe = minArray(elements.coordinates.x, numberOfElements);
	double maxXprobe = maxArray(elements.coordinates.x, numberOfElements);

	// For loop that iterate number of law that we want 
	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
        ai[iLaw] = asin((coupling.velocity / material.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI));
        ar[iLaw] = targets.tilts[iLaw] / 180 * M_PI;
        
		double a = M_PI - ai[iLaw];
		double ad = a + asin(sin(a) * (barRadius / (coupling.height + barRadius)));

		double x0 = cos((M_PI - ad)) * barRadius;
		double y1 = sin((M_PI - ad)) * barRadius;	// Y coordinate of the interface
		double x1 = (barRadius - x0);				// X coordinate of the interface

		// X coordinate of the defect
		double x2 = sin((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barRadius);
		// Y coordinate of the defect
		double y2 = cos((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barRadius);
		// Calculation of 
		double distance = sqrt(pow(x2, 2.0) + pow(y2, 2.0) + pow(0, 2.0));

		if (x0 == barRadius)
		{
			x[iLaw] = targets.positions[iLaw] + (((maxXprobe - minXprobe) / 2) + minXprobe);
			y[iLaw] = 0;

			xInt[iLaw] = 0;
			yInt[iLaw] = 0;
		}
		else
		{
			x[iLaw] = (targets.positions[iLaw] * (x2 / distance)) + (barRadius - x0) + (((maxXprobe - minXprobe) / 2) + minXprobe);;
			y[iLaw] = (targets.positions[iLaw] * (y2 / distance)) + y1;

			xInt[iLaw] = x1;
			yInt[iLaw] = y1;
		}

		z[iLaw] = targets.positions[iLaw] * (0 / distance);

		zInt[iLaw] = 0;

		utAngle[iLaw] = ar[iLaw] / M_PI * 180;
	}

	free(ai);
	free(ar);

	std::vector<double*> values{x, y, z, utAngle, xInt, yInt, zInt};

	return values;
}

std::vector<double*> Cbar::notcheBuilder(double barDiameter2)
{
	double* ai = (double*)malloc(numberOfTargets * sizeof(double));
	double* ar = (double*)malloc(numberOfTargets * sizeof(double));
	double* utAngle = (double*)malloc(numberOfTargets * sizeof(double));

	double* x = (double*)malloc(numberOfTargets * sizeof(double));
	double* y = (double*)malloc(numberOfTargets * sizeof(double));
	double* z = (double*)malloc(numberOfTargets * sizeof(double));

	double* xInt = (double*)malloc(numberOfTargets * sizeof(double));
	double* yInt = (double*)malloc(numberOfTargets * sizeof(double));
	double* zInt = (double*)malloc(numberOfTargets * sizeof(double));

	double minXprobe = minArray(elements.coordinates.x, numberOfElements);
	double maxXprobe = maxArray(elements.coordinates.x, numberOfElements);


	double* focalLength = (double*)malloc(numberOfTargets * sizeof(double));
	
	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
		if (angleType == ANGLE_TYPE::TRANSMITED)
		{
			ai[iLaw] = asin((coupling.velocity / material.velocity) * sin(targets.notchesAngles[iLaw] / 180 * M_PI));
			ar[iLaw] = targets.notchesAngles[iLaw] / 180 * M_PI;
		}

		else if (angleType == ANGLE_TYPE::INCIDENT)
		{
			ar[iLaw] = asin((material.velocity / coupling.velocity) * sin(targets.notchesAngles[iLaw] / 180 * M_PI));
			ai[iLaw] = targets.notchesAngles[iLaw] / 180 * M_PI;
		}
		
		

		double a = M_PI - ai[iLaw];
		double ad = asin(sin(a) * (barDiameter2 / (coupling.height + barDiameter2)));

		double b = M_PI - (ad + a);
		double y1 = sin(b) * barDiameter2;
		double x0 = cos(b) * barDiameter2;
		double x1 = barDiameter2 - x0;
		double d = M_PI - (((M_PI / 2) - b) + ar[iLaw]);
		double t = M_PI - (2 * ar[iLaw]);
		double h = (sin(t) / sin(ar[iLaw])) * barDiameter2;

		focalLength[iLaw] = h;
        utAngle[iLaw] = ar[iLaw] / M_PI * 180;
		
		

		double x2 = h * sin(d);
		double y2 = h * cos(d);

		x[iLaw] = x2 + x1 + (((maxXprobe - minXprobe) / 2) + minXprobe);
		y[iLaw] = y2 + y1;
		z[iLaw] = 0;

		xInt[iLaw] = x1;
		yInt[iLaw] = y1;
		zInt[iLaw] = 0;
	}

	free(ai);
	free(ar);

	std::vector<double*> values{x, y, z, utAngle, focalLength, xInt, yInt, zInt};

	return values;
}

int Cbar::Calculate()
{
	double maxXprobe = maxArray(elements.coordinates.x, numberOfElements);
	double minXprobe = minArray(elements.coordinates.x, numberOfElements);
	double maxZprobe = maxArray(elements.coordinates.z, numberOfElements);
	double minZprobe = minArray(elements.coordinates.z, numberOfElements);

	for (int iElem = 0; iElem < numberOfElements; iElem++)
	{
		elements.coordinates.y[iElem] = coupling.height - elements.coordinates.y[iElem];
	}

	// This if case handle when the defect type selected is FBH.
    if (defectType == DEFECT_TYPE::FBH){
        double* incidentAngle = (double*)malloc(numberOfTargets * sizeof(double));
        double* zDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* xDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* yDef = (double*)malloc(numberOfTargets * sizeof(double));


        for (size_t i = 0; i < numberOfTargets; i++)
        {
            incidentAngle[i] = asin(sin(targets.tilts[i] / 180 * M_PI) * (coupling.velocity / material.velocity));
        }

        double maxAngle = maxArray(incidentAngle, numberOfTargets);
        double minAngle = minArray(incidentAngle, numberOfTargets);

		free(incidentAngle);

        std::vector<double*> fbhValues = fbhBuilder(barDiameter/2);
        
        for (int i = 0; i < numberOfTargets; i++)
        {
            zDef[i] = fbhValues[2][i];
            yDef[i] = fbhValues[0][i] + coupling.height;
            xDef[i] = fbhValues[1][i];
        }

        double if1;
        double if2;

        if (tan(maxAngle) * coupling.height == tan(minAngle) * coupling.height)
        {
            if (tan(maxAngle) * coupling.height >= 0)
            {
                if1 = tan(maxAngle) * coupling.height;
                if2 = 0;
            }
            else
            {
                if1 = 0;
                if2 = tan(minAngle) * coupling.height;
            }
        }
        else
        {
            if1 = tan(maxAngle) * coupling.height;
            if2 = tan(minAngle) * coupling.height;
        }

        std::vector<double> preXIntB;
        std::vector<double> preYIntB;

		// Find the min and max element of X and Y arrays
        double minXProbe = minArray(elements.coordinates.x, numberOfElements);
        double maxXProbe = maxArray(elements.coordinates.x, numberOfElements);
        double minZProbe = minArray(elements.coordinates.z, numberOfElements);
        double maxZProbe = maxArray(elements.coordinates.z, numberOfElements);

        if (numberOfElements <= 1)
            if1 = resolution;
        if (numberOfElements > 1)
            if2 = abs(if2);
        else
            if2 = resolution;

		// This for loop is here to compute the interface points on X and Y axis
        for (int i = 0; i < ((barDiameter * M_PI / 2) / resolution) + 1; i++)
        {
            if (cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                >= minXProbe - if2 &&
                cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                < maxXProbe + if1)
            {
                preYIntB.push_back(sin(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2) + (barDiameter / 2) + coupling.height);

                preXIntB.push_back(cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2));
            }
        }

		std::vector<double> xIntB;
		std::vector<double> yIntB;
        double* zIntB = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
        * sizeof(double));

		// This for loop is here to duplicate the interface points coordinates on X and Y axis and to calculate
		// Z coordinate.
        for (int i = 0; i < (((maxZProbe - minZProbe) / resolution) + 1); i++)
        {
			for (int j = 0; j < preXIntB.size(); j++)
			{
				xIntB.push_back(preXIntB[j]);
				yIntB.push_back(preYIntB[j]);
			}
            

            for (int j = i * (preYIntB.size()); j < preYIntB.size() * (i + 1); j++)
            {
                zIntB[j] = (resolution * i) + minZProbe;
            }
        }
		
		// If _OPTIMIZATION is defined in the header file (Cbar.h), we are going to enter in this case
		// where we are going to find the minimum time for he US to go from the point to the focal point.
		// This optimized version uses an modified version of dichotomy that work in parabolas.
		#ifdef _OPTIMIZATION
		// This variable is the number of parabola that we will have.
		int nbGroupInt = (((maxZProbe - minZProbe) / resolution) + 1);

        for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
			double *compar = (double *)malloc(numberOfElements * sizeof(double));

			int decalage = xIntB.size() / nbGroupInt;

			double addTimeElemIntDef = 0;


            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
            * sizeof(double));

            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (material.velocity / 1000);
            }

			for (int i = 0; i < numberOfElements; i++)
			{
				compar[i] = INFINITY;
			}

			// changement de la boucle ici
            for (int i = 0; i < nbGroupInt; i++)
            {
				for (int probeElem = 0; probeElem < numberOfElements; probeElem++)
				{
					double* distIntProbe = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
					* sizeof(double));

					distIntProbe[0] = sqrt(pow(elements.coordinates.x[i] - xIntB[0], 2.0) 
						+ pow(elements.coordinates.y[i] - yIntB[0], 2.0) 
						+ pow(elements.coordinates.z[i] - zIntB[0], 2.0)) / (coupling.velocity / 1000);


					// Indexes and length for the while loop below.
					int start = (decalage * i);
					int end = (decalage * (i + 1)) - 1;
					int mid = (int)((end - start) / 2);
					int length = 0;

					while (length != -1)
					{
						length = end - start;
						mid = start + (int)(length / 2);

						if (mid == 0 || mid == xIntB.size())
							break;

						// Check if the previous time taken by the US is lower than the current.
						if (distDefInt[mid - 1] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[mid - 1], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[mid - 1], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[mid - 1], 2.0))) / (coupling.velocity / 1000) <
							distDefInt[mid] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[mid], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[mid], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[mid], 2.0))) / (coupling.velocity / 1000))
						{
							end = mid - 1;
						}

						// If the previous time taken by the US is bigger than the current we enter in this case
						else
						{
							start = mid + 1;
						}
					}

					addTimeElemIntDef = distDefInt[end] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[end], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[end], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[end], 2.0))) / (coupling.velocity / 1000);

					if (compar[probeElem] > addTimeElemIntDef)
					{
						compar[probeElem] = addTimeElemIntDef;
					}

					free(distIntProbe);
				}
            }
			// Maximum element of delayLaw array.
			double maxDelayLaw = maxArray(compar, numberOfElements);
			// For loop that push the delay law for each element of the probe
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				if (maxDelayLaw - compar[iElem] == NAN)
				{
					return PLUGIN_INVALID_ANGLE;
				}
				laws[iLaw].delays[iElem] = maxDelayLaw - compar[iElem];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = (((maxXprobe - minXprobe) / 2) + minXprobe);
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2] - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = (((maxZprobe - minZprobe) / 2) + minZprobe);
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = (((maxXprobe - minXprobe) / 2) + minXprobe);
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = ((elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2) - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = (((maxZprobe - minZprobe) / 2) + minZprobe);
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Release of the memory taken by delayLaw array.
			free(compar);

			free(distDefInt);
        }
		#endif

		// If _OPTIMIZATION is not defined in the header file we are going to use the non-optimized version of the algorithm
		// to find the minimum time for the US to go from the probe element to the focal
		#ifndef _OPTIMIZATION
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
            * sizeof(double));


            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (material.velocity / 1000);
            }

            double* minElems = (double*)malloc(numberOfElements * sizeof(double));

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                double* distIntProbe = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
                * sizeof(double));
                double* addDist = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
                * sizeof(double));


                for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
                {
                    distIntProbe[iIntPoint] = sqrt(pow(elements.coordinates.x[iElem] - xIntB[iIntPoint], 2.0) 
                    + pow(elements.coordinates.y[iElem] - yIntB[iIntPoint], 2.0) 
                    + pow(elements.coordinates.z[iElem] - zIntB[iIntPoint], 2.0)) / (coupling.velocity / 1000);

                    addDist[iIntPoint] = distIntProbe[iIntPoint] + distDefInt[iIntPoint];
                }

                minElems[iElem] = minArray(addDist, (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size());
            }

            double maxMinElems = maxArray(minElems, numberOfElements);

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                laws[iLaw].delays[iElem] = maxMinElems - minElems[iElem];
            }

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = elements.coordinates.x[numberOfElements / 2];
				paths[iLaw].x[1] = fbhValues[4][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2];
				paths[iLaw].y[1] = fbhValues[5][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = elements.coordinates.z[numberOfElements / 2];
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = (elements.coordinates.x[numberOfElements / 2] + elements.coordinates.x[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].x[1] = fbhValues[4][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = (elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].y[1] = fbhValues[5][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = (elements.coordinates.z[numberOfElements / 2] + elements.coordinates.z[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

            free(minElems);
        }
		#endif

		// Release of the memory taken by x/y/zDef and zIntB
		free(xDef);
		free(yDef);
		free(zDef);
		
		free(zIntB);
		
        for (size_t i = 0; i < fbhValues.size(); i++)
		{
			free(fbhValues[i]);
		}
    }

	// This else handle the case when the defect type selected is Notche.
    else
    {
        double* asinNotcheRad = (double*)malloc(numberOfTargets * sizeof(double));
        double* zDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* xDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* yDef = (double*)malloc(numberOfTargets * sizeof(double));

        for (int i = 0; i < numberOfTargets; i++)
        {
            if (angleType == ANGLE_TYPE::TRANSMITED)
            {
                asinNotcheRad[i] = asin(sin(targets.notchesAngles[i] / 180 * M_PI) * (coupling.velocity / material.velocity));
            }
            else
            {
                asinNotcheRad[i] = targets.notchesAngles[i] / 180 * M_PI;
            }
        }
		
		// Find the minimum and maximum elements of asinNotcheRad array
        double maxAngle = maxArray(asinNotcheRad, numberOfTargets);
        double minAngle = minArray(asinNotcheRad, numberOfTargets);
		
		free(asinNotcheRad);
        std::vector<double*> notcheValues = notcheBuilder(barDiameter/2);

		// For loop where we get x,y,z coordinates of notches
        for (int i = 0; i < numberOfTargets; i++)
        {
            zDef[i] = notcheValues[2][i];
            yDef[i] = notcheValues[0][i] + coupling.height;
            xDef[i] = notcheValues[1][i];
        }

        double* deflexionAngle = (double*)malloc(numberOfTargets * sizeof(double));

		// For loop that calculate the deflexion angle following the angle type
        for (int i = 0; i < numberOfTargets; i++)
        {
            if (angleType == ANGLE_TYPE::TRANSMITED)
            {
                deflexionAngle[i] = asin((sin(targets.notchesAngles[i] / 180 * M_PI) / material.velocity) * coupling.velocity);
            }
            else
            {
                deflexionAngle[i] = targets.notchesAngles[i] / 180 * M_PI;
            }

            deflexionAngle[i] = round(asin(sin(M_PI - deflexionAngle[i]) * ((barDiameter / 2) / ((barDiameter / 2) + coupling.height)))
                                / M_PI * 180);
        }

        
        double if1;
        double if2;

        if (tan(maxAngle) * coupling.height == tan(minAngle) * coupling.height)
        {
            if (tan(maxAngle) * coupling.height >= 0)
            {
                if1 = tan(maxAngle) * coupling.height;
                if2 = 0;
            }
            else
            {
                if1 = 0;
                if2 = tan(minAngle) * coupling.height;
            }
        }
        else
        {
            if1 = tan(maxAngle) * coupling.height;
            if2 = tan(minAngle) * coupling.height;
        }

        std::vector<double> preXIntB;
        std::vector<double> preYIntB;
		//inv ici
        double minXProbe = minArray(elements.coordinates.x, numberOfElements);
        double maxXProbe = maxArray(elements.coordinates.x, numberOfElements);
        double minZProbe = minArray(elements.coordinates.z, numberOfElements);
        double maxZProbe = maxArray(elements.coordinates.z, numberOfElements);

        if (numberOfElements <= 1)
            if1 = resolution;
        if (numberOfElements > 1)
            if2 = abs(if2);
        else
            if2 = resolution;


        for (int i = 0; i < ((barDiameter * M_PI / 2) / resolution) + 1; i++)
        {
            if (cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                >= minXProbe - if2 &&
                cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                < maxXProbe + if1)
            {
				//inv ici
                preYIntB.push_back(sin(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2) + (barDiameter / 2) + coupling.height);

                preXIntB.push_back(cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2));
            }
        }

		std::vector<double> xIntB;
		std::vector<double> yIntB;
		double* zIntB = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
		* sizeof(double));


        for (int i = 0; i < (((maxZProbe - minZProbe) / resolution) + 1); i++)
        {
			for (int j = 0; j < preXIntB.size(); j++)
			{
				xIntB.push_back(preXIntB[j]);
				yIntB.push_back(preYIntB[j]);
			}

            for (int j = i * (preYIntB.size()); j < preYIntB.size() * (i + 1); j++)
            {
                zIntB[j] = (resolution * i) + minZProbe;
            }
        }

		#ifdef _OPTIMIZATION
		int nbGroupInt = (((maxZProbe - minZProbe) / resolution) + 1);
		

        for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
			double *compar = (double *)malloc(numberOfElements * sizeof(double));

			// changer decalage données brut
			int decalage = xIntB.size() / nbGroupInt;

			double addTimeElemIntDef = 0;


            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1 + 1) * preYIntB.size() * sizeof(double));

            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (material.velocity / 1000);
            }


			for (int i = 0; i < numberOfElements; i++)
			{
				compar[i] = INFINITY;
			}

			// changement de la boucle ici
            for (int i = 0; i < nbGroupInt; i++)
            {
				for (int probeElem = 0; probeElem < numberOfElements; probeElem++)
				{
					double* distIntProbe = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1 + 1) * preYIntB.size()
					* sizeof(double));

					distIntProbe[0] = sqrt(pow(elements.coordinates.x[i] - xIntB[0], 2.0) 
						+ pow(elements.coordinates.y[i] - yIntB[0], 2.0) 
						+ pow(elements.coordinates.z[i] - zIntB[0], 2.0)) / (coupling.velocity / 1000);
			

					// Indexes and length for the while loop below.
					int start = (decalage * i);
					int end = (decalage * (i + 1)) - 1;
					int mid = (int)((end - start) / 2);
					int length = 0;

					while (length != -1)
					{
						length = end - start;
						mid = start + (int)(length / 2);

						// Check if the previous time taken by the US is lower than the current.
						if (distDefInt[mid - 1] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[mid - 1], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[mid - 1], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[mid - 1], 2.0))) / (coupling.velocity / 1000) <
							distDefInt[mid] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[mid], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[mid], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[mid], 2.0))) / (coupling.velocity / 1000))
						{
							end = mid - 1;
						}

						// If the previous time taken by the US is bigger than the current we enter in this case
						else
						{
							start = mid + 1;
						}
					}


					addTimeElemIntDef = distDefInt[end] + (sqrt(pow(elements.coordinates.x[probeElem] - xIntB[end], 2.0) + pow(elements.coordinates.y[probeElem] - yIntB[end], 2.0) + pow(elements.coordinates.z[probeElem] - zIntB[end], 2.0))) / (coupling.velocity / 1000);

					if (compar[probeElem] > addTimeElemIntDef)
						{
							compar[probeElem] = addTimeElemIntDef;
						}
					
					free(distIntProbe);

				}
            }
			
			// Maximum element of delayLaw array.
			double maxDelayLaw = maxArray(compar, numberOfElements);
			// For loop that push the delay law for each element of the probe
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				laws[iLaw].delays[iElem] = maxDelayLaw - compar[iElem];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
				if (numberOfElements % 2 != 0) {
					paths[iLaw].x[0] = elements.coordinates.x[numberOfElements / 2];
					paths[iLaw].x[1] = notcheValues[6][iLaw];
					paths[iLaw].x[2] = xDef[iLaw];

					paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2] - coupling.height;
					paths[iLaw].y[1] = notcheValues[5][iLaw];
					paths[iLaw].y[2] = yDef[iLaw];

					paths[iLaw].z[0] = elements.coordinates.z[numberOfElements / 2];
					paths[iLaw].z[1] = notcheValues[7][iLaw];
					paths[iLaw].z[2] = zDef[iLaw];
				}

				// Get the value of remarkable x,y,z coordinates when the number of element is peer.
				else
				{
					paths[iLaw].x[0] = (elements.coordinates.x[numberOfElements / 2] + elements.coordinates.x[(numberOfElements / 2) - 1]) / 2;
					paths[iLaw].x[1] = notcheValues[6][iLaw];
					paths[iLaw].x[2] = xDef[iLaw];

					paths[iLaw].y[0] = ((elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2) - coupling.height;
					paths[iLaw].y[1] = notcheValues[5][iLaw];
					paths[iLaw].y[2] = yDef[iLaw];

					paths[iLaw].z[0] = (elements.coordinates.z[numberOfElements / 2] + elements.coordinates.z[(numberOfElements / 2) - 1]) / 2;
					paths[iLaw].z[1] = notcheValues[7][iLaw];
					paths[iLaw].z[2] = zDef[iLaw];
				}



			// Release of the memory taken by delayLaw array.
			free(compar);
			free(distDefInt);
		}
		#endif

		#ifndef _OPTIMIZATION
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
            * sizeof(double));

            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (material.velocity / 1000);
            }

            double* minElems = (double*)malloc(numberOfElements * sizeof(double));

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                double* distIntProbe = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
                * sizeof(double));
                double* addDist = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
                * sizeof(double));


                for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
                {
                    distIntProbe[iIntPoint] = sqrt(pow(elements.coordinates.x[iElem] - xIntB[iIntPoint], 2.0) 
                    + pow(elements.coordinates.y[iElem] - yIntB[iIntPoint], 2.0) 
                    + pow(elements.coordinates.z[iElem] - zIntB[iIntPoint], 2.0)) / (coupling.velocity / 1000);

                    addDist[iIntPoint] = distIntProbe[iIntPoint] + distDefInt[iIntPoint];
                }

                minElems[iElem] = minArray(addDist, (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size());
            }

            double maxMinElems = maxArray(minElems, numberOfElements);

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                laws[iLaw].delays[iElem] = maxMinElems - minElems[iElem];
            }

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = elements.coordinates.x[numberOfElements / 2];
				paths[iLaw].x[1] = notcheValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2];
				paths[iLaw].y[1] = notcheValues[6][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = elements.coordinates.z[numberOfElements / 2];
				paths[iLaw].z[1] = notcheValues[7][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = (elements.coordinates.x[numberOfElements / 2] + elements.coordinates.x[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].x[1] = notcheValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = (elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].y[1] = notcheValues[6][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = (elements.coordinates.z[numberOfElements / 2] + elements.coordinates.z[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].z[1] = notcheValues[7][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

            free(minElems);
        }
		#endif
		
		free(xDef);
		free(yDef);
		free(zDef);

		free(deflexionAngle);

		free(zIntB);

		for (int i =0; i < notcheValues.size(); i++)
		{
			free(notcheValues[i]);
		}
    }

	// Retreive the original coordinates.
	for (int iElem = 0; iElem < numberOfElements; iElem++)
	{
		elements.coordinates.y[iElem] = -1.0 * elements.coordinates.y[iElem] + coupling.height;
	}


	return PLUGIN_NO_ERROR;
}
