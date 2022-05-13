#include "framework.h"
#include "Cbar.h"
#include "unit.h"
#include "error_codes.h"
#include <cstddef>
#include <string.h>
#include "framework.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#define _CRT_SECURE_NO_WARNINGS

#define M_PI 3.14159265358979323846264338327950288

#define strcmpi(x, y) strcasecmp((x), (y))
#define _strnicmp(x, y, z) strncasecmp((x), (y), (z))

// Ajouter flag /fp:fast et /O2 a la compilation.

// TODO : ajouter les getter / setter

Cbar::Cbar()
{
	// Initialize Cplate's members

	opened = false;
	calculationDone = false;

	resolution = 0.1;
	double barDiameter;
    double* barNotchAngle;

	numberOfElements = 0;

	elements.coordinates.x = NULL;
	elements.coordinates.y = NULL;
	elements.coordinates.z = NULL;

	material.velocity = 3900.0;

	coupling.velocity = 1500.0;
	coupling.height = 0.0;

	angleType = ANGLE_TYPE::INCIDENT;

	numberOfTargets = 0;

	targets.tilts = NULL;
	targets.positions = NULL;
	laws = NULL;
	paths = NULL;

	utAngle = 1;
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

	if (strcmpi(param_name, "Flow.PositionType") == 0)
	{
		angleType = (ANGLE_TYPE)*value;
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

	if (strcmpi(param_name, "Flow.PositionType") == 0)
	{
		angleType = (ANGLE_TYPE)*value;
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
			targets.tilts[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_rad);
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

int Cbar::minArrayIndex(const double *array, int size)
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

double *Cbar::append(double *ar1, double *ar2, int len1, int len2)
{
	double *arTmp = (double *)malloc(sizeof(double) * (len1 + len2));

	for (int i = 0; i < len1; i++)
	{
		arTmp[i] = ar1[i];
	}

	int cpt = 0;

	for (int i = len1; i < len1 + len2; i++)
	{
		arTmp[i] = ar2[cpt];
		cpt++;
	}

	free(ar1);

	return arTmp;
}

std::vector<double*> Cbar::fbhBuilder(double* tiltRad, double barDiameter2)
{
	double* ai = (double*)malloc(numberOfTargets * sizeof(double));
	double* ar = (double*)malloc(numberOfTargets * sizeof(double));

	double* utAngle = (double*)malloc(numberOfTargets * sizeof(double));

	double* x = (double*)malloc(numberOfTargets * sizeof(double));
	double* y = (double*)malloc(numberOfTargets * sizeof(double));
	double* z = (double*)malloc(numberOfTargets * sizeof(double));

	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
		if (angleType == ANGLE_TYPE::TRANSMITED)
		{
			ai[iLaw] = asin((material.velocity / coupling.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI));
			ar[iLaw] = targets.tilts[iLaw] / 180 * M_PI;
		}

		else if (angleType == ANGLE_TYPE::INCIDENT)
		{
			ar[iLaw] = asin((coupling.velocity / material.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI));
			ai[iLaw] = targets.tilts[iLaw] / 180 * M_PI;
		}

		double a = ai[iLaw] - M_PI;
		double ad = a + (sin(a) * (barDiameter2 / (coupling.height + barDiameter2)));
		double y1 = sin((M_PI - ad)) * barDiameter2;
		double x0 = cos((M_PI - ad)) * barDiameter2;
		double x1 = (barDiameter2 - x0);

		double x2 = sin((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barDiameter);
		double y2 = cos((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barDiameter);

		double distance = sqrt(pow(x2, 2.0) + pow(y2, 2.0) + pow(0, 2.0));

		if (x0 == barDiameter2)
		{
			x[iLaw] = targets.positions[iLaw];
			y[iLaw] = 0;
		}
		else
		{
			x[iLaw] = (targets.positions[iLaw] * (x2 / distance)) + (barDiameter2 - x0);
			y[iLaw] = (targets.positions[iLaw] * (y2 / distance)) + y1;
		}

		z[iLaw] = targets.positions[iLaw] * (0 / distance);


		utAngle[iLaw] = ar[iLaw] / M_PI * 180;
		
		
	}

	std::vector<double*> values{x, y, z, utAngle};

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

		double x2 = h * sin(d);
		double y2 = h * cos(d);

		x[iLaw] = x2 + x1;
		y[iLaw] = y2 + y1;
		z[iLaw] = 0;
	}

	std::vector<double*> values{x, y, z, utAngle, focalLength};

	return values;
}

// nbr de notches = numberoftargets




int Cbar::Calculate()
{
	


	return PLUGIN_NO_ERROR;
}
