#include "framework.h"
#include "Ctube.h"
#include "unit.h"
#include "error_codes.h"
#include "framework.h"

#define _OPTIMIZATION

#define _CRT_SECURE_NO_WARNINGS

#define strcmpi(x, y) strcasecmp((x), (y))
#define _strnicmp(x, y, z) strncasecmp((x), (y), (z))

// Ajouter flag /fp:fast et /O2 a la compilation.

Ctube::Ctube()
{
	// Initialize Ctube's members

	opened = false;
	calculationDone = false;

	numberOfElements = 0;

	diameter = 0;
	thickness = 0;
	resolution = 0.1;
	tubeOffset = 0;

	defectType = DEFECT_TYPE::NOTCHE;

	xi3D = NULL;
	yi3D = NULL;
	zi3D = NULL;
	alphaS = NULL;

	elements.coordinates.x = NULL;
	elements.coordinates.y = NULL;
	elements.coordinates.z = NULL;

	material.velocity = 3230.0;

	coupling.velocity = 1500.0;
	coupling.height = 0;

	focal.length.coupling = 0;

	numberOfTargets = 0;

	targets.tilts = NULL;
	targets.skews = NULL;
	targets.positions = NULL;
	laws = NULL;
	paths = NULL;
}

Ctube::~Ctube()
{
	Close();
}

int Ctube::Open()
{
	if (opened)
		return PLUGIN_ALREADY_OPEN;

	opened = true;
	return PLUGIN_NO_ERROR;
}

int Ctube::Close()
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
	if (alphaS != NULL)
		free(alphaS);
	if (xi3D != NULL)
		free(xi3D);
	if (yi3D != NULL)
		free(yi3D);
	if (zi3D != NULL)
		free(zi3D);
	

	return PLUGIN_NO_ERROR;
}

int Ctube::Set(const char *param_name, int unit, int *value)
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (strcmpi(param_name, "Defect_Type") == 0)
	{
		defectType = (DEFECT_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}

int Ctube::Set(const char *param_name, int unit, double *value)
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

	if (strcmpi(param_name, "Focal.Length.Coupling") == 0)
	{
		focal.length.coupling = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Thickness") == 0)
	{
		thickness = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Diameter") == 0)
	{
		diameter = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Defect_Type") == 0)
	{
		defectType = (DEFECT_TYPE)*value;
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "tubeOffset") == 0)
	{
		tubeOffset = Unit::ChangeUnit(*value, unit, UNIT_mm);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Ctube::Set(const char *param_name, int unit, char *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Ctube::Set(const char *param_name, int unit, int *dim1, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Ctube::Set(const char *param_name, int unit, int *dim1, double *value)
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
			free(targets.skews);
			free(targets.positions);
			free(laws);
			free(paths);
			targets.tilts = NULL;
			targets.skews = NULL;
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

	if (strcmpi(param_name, "Targets.Skews") == 0)
	{
		if (targets.skews == NULL)
			targets.skews = (double *)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.skews[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_deg);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "Targets.Positions") == 0)
	{
		if (targets.positions == NULL)
			targets.positions = (double *)malloc(*dim1 * sizeof(double));
		for (int iTarget = 0; iTarget < *dim1; iTarget++)
			targets.positions[iTarget] = Unit::ChangeUnit(value[iTarget], unit, UNIT_deg);
		calculationDone = false;
		return PLUGIN_NO_ERROR;
	}

	calculationDone = false;
	return PLUGIN_UNKNOWN_PARAMETER;
}
int Ctube::Set(const char *param_name, int unit, int *dim1, int *dim2, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Set(const char *param_name, int unit, int *dim1, int *dim2, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Set(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Set(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Ctube::Get(const char *param_name, int unit, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Get(const char *param_name, int unit, double *value)
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
int Ctube::Get(const char *param_name, int unit, int *dim1, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Get(const char *param_name, int unit, int *dim1, double *value) 
{
	if (!opened)
		return PLUGIN_NOT_OPEN;

	if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;

	if (strcmpi(param_name, "AlphaS") == 0)
	{

		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			value[iLaw] = Unit::ChangeUnit(alphaS[iLaw], UNIT_deg, unit);
		}
		return PLUGIN_NO_ERROR;
	}
		
	if (strcmpi(param_name, "xi3D") == 0)
	{
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			value[iLaw] = Unit::ChangeUnit(xi3D[iLaw], UNIT_mm, unit);
		}
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "yi3D") == 0)
	{
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			value[iLaw] = Unit::ChangeUnit(yi3D[iLaw], UNIT_mm, unit);
		}
		return PLUGIN_NO_ERROR;
	}

	if (strcmpi(param_name, "zi3D") == 0)
	{
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			value[iLaw] = Unit::ChangeUnit(zi3D[iLaw], UNIT_mm, unit);
		}
		return PLUGIN_NO_ERROR;
	}


	return PLUGIN_UNKNOWN_PARAMETER; 
}
int Ctube::Get(const char *param_name, int unit, int *dim1, char *value)
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
int Ctube::Get(const char *param_name, int unit, int *dim1, int *dim2, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Get(const char *param_name, int unit, int *dim1, int *dim2, double *value)
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
		
		int nbPointsRemarkable = 4;

		if (defectType == DEFECT_TYPE::FBH)
			nbPointsRemarkable = 3;
		

		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < nbPointsRemarkable)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < nbPointsRemarkable; iPoint++)
			{
				*(value + iTarget * nbPointsRemarkable + iPoint) = Unit::ChangeUnit(paths[iTarget].x[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}
	if (strcmpi(param_name, "Paths.y") == 0)
	{
		int nbPointsRemarkable = 4;

		if (defectType == DEFECT_TYPE::FBH)
			nbPointsRemarkable = 3;

		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < nbPointsRemarkable)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < nbPointsRemarkable; iPoint++)
			{
				*(value + iTarget * nbPointsRemarkable + iPoint) = Unit::ChangeUnit(paths[iTarget].y[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}
	if (strcmpi(param_name, "Paths.z") == 0)
	{
		int nbPointsRemarkable = 4;

		if (defectType == DEFECT_TYPE::FBH)
			nbPointsRemarkable = 3;

		if (!calculationDone)
			return PLUGIN_DATA_NOT_READY;
		if (*dim1 < numberOfTargets)
			return PLUGIN_DIM1_TOO_SMALL;
		if (*dim2 < nbPointsRemarkable)
			return PLUGIN_DIM2_TOO_SMALL;

		for (int iTarget = 0; iTarget < numberOfTargets; iTarget++)
		{
			for (int iPoint = 0; iPoint < nbPointsRemarkable; iPoint++)
			{
				*(value + iTarget * nbPointsRemarkable + iPoint) = Unit::ChangeUnit(paths[iTarget].z[iPoint], UNIT_mm, unit);
			}
		}
		return PLUGIN_NO_ERROR;
	}

	return PLUGIN_UNKNOWN_PARAMETER;
}
int Ctube::Get(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, int *value) { return PLUGIN_UNKNOWN_PARAMETER; }
int Ctube::Get(const char *param_name, int unit, int *dim1, int *dim2, int *dim3, double *value) { return PLUGIN_UNKNOWN_PARAMETER; }

int Ctube::ExecSync(const char *action)
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
		if (numberOfTargets == 0 || targets.tilts == NULL || targets.skews == NULL || diameter == 0 || thickness == 0 || focal.length.coupling == 0)
			return PLUGIN_PARMETER_UNDEFINED;

		// Allocate resources for laws and paths
		laws = (Ctube::PLAW)malloc(numberOfTargets * sizeof(Ctube::LAW));
		paths = (Ctube::PPATH)malloc(numberOfTargets * sizeof(Ctube::PATH));
		
		xi3D = (double*)malloc(numberOfTargets * sizeof(double));
		yi3D = (double*)malloc(numberOfTargets * sizeof(double));
		zi3D = (double*)malloc(numberOfTargets * sizeof(double));
		alphaS = (double*)malloc(numberOfTargets * sizeof(double));

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

int Ctube::ExecAsync(const char *action) { return PLUGIN_UNKNOWN_PARAMETER; }

/**
 * maxArray
 *
 * get the max element from a array
 *
 * @param array<double> array : The array you want to get the max element
 *
 * @return double max : The max element of the array
 **/
double Ctube::maxArray(const double *array, int size)
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
 * @param size : the size of the array
 *
 * @return double min : The min element of the array
 **/
double Ctube::minArray(const double *array, int size)
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



/**
 * @brief This fonction return all remarkable point of the US in the tube and the deflection angle.
 * 
 * @param skew 
 * @param alphaI 
 * @return std::vector<double> 
 */
std::vector<double> Ctube::newElipse(double skew, double alphaI)
{	
	using namespace std;		// Using namespace here to simplify the code.
	using namespace Eigen;

	double penetrationAngle;
	double alphaR = 90.0 - asin(((material.velocity / coupling.velocity) * sin(alphaI / 180.0 * M_PI))) / M_PI * 180.0;

	penetrationAngle = alphaR;
	
	// sinus and cosinus of gamma angle.
	double cga = cos(-M_PI/2.0);
	double sga = sin(-M_PI/2.0);	

	MatrixXf mat(4,4);
	mat << cga, -sga, 0, 0 ,
	sga, cga, 0, 0,
	0, 0 ,1 ,0,
	0, 0, 0, 1;


	double alphaObl = skew / 180.0 * M_PI;
	double betaProf = -((penetrationAngle / 180.0) * M_PI);

	// Calculation of the director vector M1_M2
	MatrixXf dirD1D2(4,1);
	dirD1D2 << cos(alphaObl) * cos(betaProf), sin(betaProf), cos(betaProf) * sin(alphaObl), 0.0;

	MatrixXf dir2;
	dir2 = mat * dirD1D2;

	// Direction of the default in the calculation landmark
	double teta = (atan2(dir2(1), dir2(0)) * 180.0) / M_PI;
	double alpha = (atan2(dir2(2), dir2(0)) * 180.0) / M_PI;

	
	double r1 = diameter / 2.0;		// Radius of the inner wall.
	double r2 = r1 - thickness;		// Radius of the outer wall.

	// x,y,z coordinates  of the defect
	double xm1 = r1;
	double ym1 = 0.0;
	double zm1 = 0.0;

	// Matrix of coordinate of the defect.
	MatrixXf m1(4,1);
	m1 << r1, 0.0, 0.0, 1.0;

	alpha = alpha / 180.0 * M_PI;
	teta = (teta / 180.0) * M_PI;

	// Equation of the segment M1_M2 (M1 is the defect and M2 is the point where the US hit the
	// inner wall).
	double beta = atan(cos(alpha) * tan(teta));
	double a1 = cos(alpha) * cos(beta);
	double b1 = sin(beta);
	double c1 = cos(beta) * std::sin(alpha);

	// Calculation of the discriminant to solve the equation above.
	double discri = sqrt(pow((xm1 * 2) * a1, 2.0) - ((pow(a1, 2.0) + pow(b1, 2.0)) * (pow(r1, 2.0) - pow(r2, 2.0)) * 4)) + -((xm1 * 2) * a1);

	// tPp and tP are the 2 different solution that we can obtain because its a quadratic equation.
	double tPp = discri / ((pow(a1, 2.0) + pow(b1, 2.0)) * 2);
	double tP = (-((xm1 * 2.0) * a1) - 
	sqrt(pow((xm1 * 2) * a1, 2.0) - ((pow(a1, 2.0) + pow(b1, 2.0)) * 
	(pow(r1, 2.0) - pow(r2, 2.0)) * 4))) / 
	((pow(a1, 2.0) + pow(b1, 2.0)) * 2);

	// Here is the variable where we are going to put the shortest path between tPp and tP.
	double T_LgM1M2;

	if (tPp < tP)
		T_LgM1M2 = tPp;
	else
		T_LgM1M2 = tP;
	
	// Here we have the x,y,z coordinates of where the US hit the inner wall of the tube.
	double xm2 = (T_LgM1M2 * a1) + xm1;
	double ym2 = b1 * T_LgM1M2;
	double zm2 = T_LgM1M2 * c1;

	// Matrix of M2 coordinates.
	MatrixXf m2(4,1);
	m2 << xm2, ym2, zm2, 1;

	double gamma2 = atan(ym2 / xm2);

	// Calculation of sinus and cosinus of gamma angle for m2.
	double sgam = sin(gamma2);
	double cgam = cos(gamma2);

	// Matrix for the rotation-translation to change to the M2 landmark.
	Matrix4f mat1;
	mat1 << cgam, sgam, 0, -r2,
	-sgam, cgam, 0, 0,
	0, 0, 1.0, -zm2,
	0, 0, 0, 1.0;

	// Matrix coordinates of M1 in the M2 landmark.
	MatrixXf m1p;
	m1p = mat1 * m1;

	// Matrix of coordinate of the reflexion point of M1.
	MatrixXf m1pr(4, 1);
	m1pr << m1p(0), -m1p(1), -m1p(2), 1;

	// Coordinates of the point reflected from M1 by the tangent plane.
	MatrixXf m1r;
	m1r = mat1.inverse() * m1pr;

	// Distance from M2 to M1R.
	double lg_m2_m1r = sqrt(pow(xm2 - m1r(0), 2.0) + pow(ym2 - m1r(1), 2.0) + pow(zm2 - m1r(2), 2.0));

	// Equation of the segment M1R to M2.
	double a2 = (xm2 - m1r(0)) / lg_m2_m1r;
	double b2 = (ym2 - m1r(1)) / lg_m2_m1r;
	double c2 = (zm2 - m1r(2)) / lg_m2_m1r;

	// Calculation of the discriminant to solve the equation above.
	double discriminant = (pow((ym2 * 2 * b2) + (xm2 * 2 * a2),2.0)) - 
	(((pow(xm2, 2.0) + pow(ym2, 2.0)) - pow(xm1, 2.0)) * (pow(a2, 2.0) + pow(b2, 2.0)) * 4);

	// The 2 solution to the equation.
	double tPp2 = (sqrt(discriminant) + -((ym2 * 2 * b2) + (xm2 * 2 * a2))) / ((pow(a2, 2.0) + pow(b2, 2.0)) * 2);
	double tP2 = (-((ym2 * 2 * b2) + (xm2 * 2 * a2)) - sqrt(discriminant)) / ((pow(a2, 2.0) + pow(b2, 2.0)) * 2);

	// Here is the variable where we are going to put the shortest path between tPp2 and tP2.
	double t_lg_m2m3;

	if (tPp2 < tP2)
		t_lg_m2m3 = tPp2;
	else
		t_lg_m2m3 = tP2;

	// x,y,z coordinates of M3 point.
	double xm3 = a2 * t_lg_m2m3 + xm2;
	double ym3 = b2 * t_lg_m2m3 + ym2;
	double zm3 = c2 * t_lg_m2m3 + zm2;

	// Matrix of M3 point coordinates.
	MatrixXf m3(4, 1);
	m3 << xm3, ym3, zm3, 1;

	double gamma3 = atan(ym3 / xm3);

	double sgam2 = sin(gamma3);
	double cgam2 = cos(gamma3);

	// Matrix for the rotation-translation to change to the M3 landmark.
	Matrix4f mat3;
	mat3 << cgam2, sgam2, 0, -xm1,
		-sgam2, cgam2, 0, 0,
		0, 0, 1, -zm3,
		0, 0, 0, 1;

	// Coordinates of M2 point in the M3 landmark.
	MatrixXf m2s;
	m2s = mat3 * m2;

	// Distance between M2 and M3 points.
	double lgM2M3 = sqrt(pow(m2s(0), 2.0) + pow(m2s(1), 2.0) + pow(m2s(2), 2.0));
	double lget = sqrt(pow(m2s(1), 2.0) + pow(m2s(2), 2.0));
	
	// Director vector M2_M3.
	MatrixXf dirM2M3(4, 1);
	dirM2M3 << -m2s(0) / lgM2M3, -m2s(1) / lgM2M3, -m2s(2) / lgM2M3, 0;

	// aglT is the angle transmited in the material.
	double aglT = (acos(-(m2s(0) / lgM2M3)) / M_PI) * 180;
	// aglI is the angle transmited in the coupling.
	double aglI = (asin(sin(acos(-(m2s(0,0) / lgM2M3))) * (coupling.velocity / material.velocity)) / M_PI * 180);

	// xI,yI,zI are the coordinates of the incident vector in the coupling.
	double xI = -m2s(0);

	double lgei = abs(tan(asin(sin(acos(-(m2s(0) / lgM2M3))) * (coupling.velocity / material.velocity))) * xI);

	double yI = (lgei / lget) * -m2s(1);
	double zI = (lgei / lget) * -m2s(2);

	// Matrix of coordinates of the incident vector in the coupling
	MatrixXf m2ts(4,1);
	m2ts << xI, yI, zI, 1;

	// Coordinates of M2t point in the 0xyz landmark.
	MatrixXf m2t;
	m2t = mat3.inverse() * m2ts;

	// The exact same method as to get M1-M2 path.
	double r4 = xm1 + coupling.height;

	double lgM3M2t = sqrt(pow(xm3 - m2t(0), 2.0) + pow(ym3 - m2t(1), 2.0) + pow(zm3 - m2t(2), 2.0));

	double a3 = (xm3 - m2t(0)) / lgM3M2t;
	double b3 = (ym3 - m2t(1)) / lgM3M2t;
	double c3 = (zm3 - m2t(2)) / lgM3M2t;

	double discriminant2 = sqrt(pow((xm3 * 2 * a3) + (ym3 * 2 * b3), 2.0) - (((pow(xm3, 2.0) + pow(ym3, 2.0)) - pow(r4, 2.0)) * (pow(a3, 2.0) + pow(b3, 2.0)) * 4));
	double tPp3 = (-((xm3 * 2 * a3) + (ym3 * 2 * b3)) + discriminant2) / ((pow(a3, 2.0) + pow(b3, 2.0)) * 2);
	double tP3 = (-((xm3 * 2 * a3) + (ym3 * 2 * b3)) - discriminant2) / ((pow(a3, 2.0) + pow(b3, 2.0)) * 2);
	
	double tLgM3M4;

	if (tPp3 < tP3)
		tLgM3M4 = tPp3;
	else
		tLgM3M4 = tP3;
	
	double xm4 = tLgM3M4 * a3 + xm3;
	double ym4 = tLgM3M4 * b3 + ym3;
	double zm4 = tLgM3M4 * c3 + zm3;

	MatrixXf m4(4, 1);
	m4 << xm4, ym4, zm4, 1;

	double gamma4 = atan(ym4 / xm4);
	double sgam4 = sin(gamma4);
	double cgam4 = cos(gamma4);

	Matrix4f mat5;
	mat5 << cgam4 , sgam4, 0, -r4,
	-sgam4, cgam4, 0, 0,
	0, 0, 1, -zm4,
	0, 0, 0, 1;

	// Coordinate x,y,z of M1,M2,M3,M4.
	MatrixXf m1PhC(4, 1);
	m1PhC = mat5 * m1;
	MatrixXf m2PhC(4, 1);
	m2PhC = mat5 * m2;
	MatrixXf m3PhC(4, 1);
	m3PhC = mat5 * m3;
	MatrixXf m4PhC(4, 1);
	m4PhC = mat5 * m4;

	// Coordinate of the defect
	// Inversion x y 
	double yDef = m1PhC(0);
	double xDef = m1PhC(1);
	double zDef = m1PhC(2);
	
	// Coordinates of the reflex in the tube
	double yRef = m2PhC(0);
	double xRef = m2PhC(1);
	double zRef = m2PhC(2);

	// Value of interface points for the paths.
	double yInter = m3PhC(0);
	double xInter = m3PhC(1);
	double zInter = m3PhC(2);

	// Absolute value of interface points.
	double yi3D = abs(m3PhC(0));
	double xi3D = abs(m3PhC(1));
	double zi3D = abs(m3PhC(2));

	// Angle of incidence on the surface.
	// inversion x y
	double alphaS = acos(yi3D / sqrt(pow(yi3D, 2.0) + pow(xi3D, 2.0) + pow(zi3D, 2.0))) / M_PI * 180;

	// Here we get all the values in a vector to use them in Calculate function.
	vector<double> returnVal {xi3D, yi3D, zi3D, alphaS, xDef, yDef, zDef, xRef, yRef, zRef, xInter, yInter, zInter};

	return returnVal;
}


std::vector<double*> Ctube::fbhBuilder(double barRadius)
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

	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
        ai[iLaw] = asin((coupling.velocity / material.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI));
        ar[iLaw] = targets.tilts[iLaw] / 180 * M_PI;
        
		double a = M_PI - ai[iLaw];
		double ad = a + asin(sin(a) * (barRadius / (coupling.height + barRadius)));

		

		double y1 = sin((M_PI - ad)) * barRadius;
		double x0 = cos((M_PI - ad)) * barRadius;
		double x1 = (barRadius - x0);

		double x2 = sin((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barRadius);
		double y2 = cos((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barRadius);

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
			x[iLaw] = (targets.positions[iLaw] * (x2 / distance)) + (barRadius - x0) + (((maxXprobe - minXprobe) / 2) + minXprobe);
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

// Calculate function return the the delay for each element of the probe for each law.
// It also return remarkable points for each law.
int Ctube::Calculate()
{	
	for (int iElem = 0; iElem < numberOfElements; iElem++)
	{
		elements.coordinates.y[iElem] = coupling.height - elements.coordinates.y[iElem];
	}

	// Case where the defect type selected is a notch.
	if (defectType == DEFECT_TYPE::NOTCHE)
	{
		double maxYProbe = maxArray(elements.coordinates.x, numberOfElements);
		double minYProbe = minArray(elements.coordinates.x, numberOfElements);

		// For loop that iterate the number of law that we want.
		for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
		{
			std::vector<double> elipse = Ctube::newElipse(targets.skews[iLaw], targets.tilts[iLaw]);

			// Here we get all the values given by the function newElipse
			double xI3Dv = elipse.at(0) + tubeOffset;
			double yI3Dv = elipse.at(1);
			double zI3Dv = elipse.at(2);
			double alphaSv = elipse.at(3);

			// XYZ coordinates of the defect.
			double xDef = elipse.at(4);
			double yDef = elipse.at(5);
			double zDef = elipse.at(6);
			// XYZ coordinates of the reflection point in the tube.
			double xRef = elipse.at(7);
			double yRef = elipse.at(8);
			double zRef = elipse.at(9);
			// XYZ coordinates of the interface point.
			double xInter = elipse.at(10);
			double yInter = elipse.at(11);
			double zInter = elipse.at(12);

			// Put x and z coordinates in the right dial by inverting x and z coordinates depending on the tilt and the skew angle asked.
			double xRightDial;
			double zRightDial;

			if (targets.skews[iLaw] <= 360 && targets.skews[iLaw] >= 180)
				zRightDial = -(zI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
			else
				zRightDial = (zI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
			
			if (targets.skews[iLaw] <= 270 && targets.skews[iLaw] >= 90)
				xRightDial = -(xI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
			else
				xRightDial = (xI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
			
			// Array of all the distances between the probe and the defect.
			double* distancesArray = (double*)malloc(numberOfElements * sizeof(double));
			
			// For loop where are going to compute all possible path.
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				double dist2 = sqrt(pow(((yI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0)))) * focal.length.coupling)
				- elements.coordinates.y[iElem], 2.0) 
				+ pow((xRightDial * focal.length.coupling) - elements.coordinates.x[iElem], 2.0)
				+ pow((zRightDial * focal.length.coupling) - elements.coordinates.z[iElem], 2.0));

				distancesArray[iElem] = dist2;
			}
					
			
			// Maximum value of tabDist array
			double maxDistancesArray = maxArray(distancesArray, numberOfElements);


			// For loop that iterate number of element times and where we are going to compute
			// all the delay for each element (converted to milliseconds).
			for (int iElem = 0; iElem < numberOfElements; iElem++)
			{
				distancesArray[iElem] = (((maxDistancesArray - distancesArray[iElem]) / coupling.velocity) * 1000);

				laws[iLaw].delays[iElem] = distancesArray[iElem];
			}
			
			// Release of the memory taken by distances array
			free(distancesArray);

			xI3Dv -= tubeOffset;

			// Here we get the values from the newElipse function for each law.
			alphaS[iLaw] = alphaSv;
			xi3D[iLaw] = xI3Dv;
			yi3D[iLaw] = yI3Dv;
			zi3D[iLaw] = zI3Dv;
			
			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = elements.coordinates.x[numberOfElements / 2];
				paths[iLaw].x[1] = xInter;
				paths[iLaw].x[2] = xRef;
				paths[iLaw].x[3] = xDef;

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2];
				paths[iLaw].y[1] = yInter;
				paths[iLaw].y[2] = yRef;
				paths[iLaw].y[3] = yDef;

				paths[iLaw].z[0] = elements.coordinates.z[numberOfElements / 2];
				paths[iLaw].z[1] = zInter;
				paths[iLaw].z[2] = zRef;
				paths[iLaw].z[3] = zDef;
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = (elements.coordinates.x[numberOfElements / 2] + elements.coordinates.x[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].x[1] = xInter;
				paths[iLaw].x[2] = xRef;
				paths[iLaw].x[3] = xDef;

				paths[iLaw].y[0] = (elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].y[1] = yInter;
				paths[iLaw].y[2] = yRef;
				paths[iLaw].y[3] = yDef;

				paths[iLaw].z[0] = (elements.coordinates.z[numberOfElements / 2] + elements.coordinates.z[(numberOfElements / 2) - 1]) / 2;
				paths[iLaw].z[1] = zInter;
				paths[iLaw].z[2] = zRef;
				paths[iLaw].z[3] = zDef;
			}

		}
	}	// End of the case where the defect type is notch.

	// Else case that handle when the defect is a FBH.
	else
	{
		double* incidentAngle = (double*)malloc(numberOfTargets * sizeof(double));
        double* zDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* xDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* yDef = (double*)malloc(numberOfTargets * sizeof(double));

		// For loop that calculate all the incident angle for each law.
        for (int i = 0; i < numberOfTargets; i++)
        {
            incidentAngle[i] = asin(sin(targets.tilts[i] / 180 * M_PI) * (coupling.velocity / material.velocity));
        }

        double maxAngle = maxArray(incidentAngle, numberOfTargets);
        double minAngle = minArray(incidentAngle, numberOfTargets);

		free(incidentAngle);

		// Here we get the values from the fbhBuilder function the get the coordinates XYZ of the fbh.
        std::vector<double*> fbhValues = fbhBuilder(diameter/2);
        
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

		// For loop where we calculate the XY coordinates of the interfaces points before duplicate it.
        for (int i = 0; i < ((diameter * M_PI / 2) / resolution) + 1; i++)
        {
            if (cos(((M_PI / ((diameter * M_PI / 2) / resolution)) * i) - M_PI) * (diameter / 2)
                >= minXProbe - if2 &&
                cos(((M_PI / ((diameter * M_PI / 2) / resolution)) * i) - M_PI) * (diameter / 2)
                < maxXProbe + if1)
            {	
                preYIntB.push_back(sin(((M_PI / ((diameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (diameter / 2) + (diameter / 2) + coupling.height);

                preXIntB.push_back(cos(((M_PI / ((diameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (diameter / 2));
            }
        }

		
		std::vector<double> xIntB;
		std::vector<double> yIntB;
        double* zIntB = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
        * sizeof(double));

		// For loop where we duplicate the XY coordinates of the interfaces points.
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

		// If _OPTIMIZATION is defined in Cplate.h, we use the optimization function to calculate the
		// total path of the US.
		#ifdef _OPTIMIZATION
		int nbGroupInt = (((maxZProbe - minZProbe) / resolution) + 1);

        for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
			double *compar = (double *)malloc(numberOfElements * sizeof(double));

			int offset = xIntB.size() / nbGroupInt;

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
					int start = (offset * i);
					int end = (offset * (i + 1)) - 1;
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
				if (maxDelayLaw - compar[iElem] == NAN)
				{
					return PLUGIN_INVALID_ANGLE;
				}
				laws[iLaw].delays[iElem] = maxDelayLaw - compar[iElem];
			}

			double deflectionAngle = (M_PI - asin((coupling.velocity / material.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI)))
			+ asin(sin(M_PI - asin((coupling.velocity / material.velocity) * sin(targets.tilts[iLaw] / 180 * M_PI)))
			* ((diameter / 2) / (coupling.height + (diameter / 2))));

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = ((maxXProbe - minXProbe) / 2) + minXProbe;
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2] - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = ((maxZProbe - minZProbe) / 2) + minZProbe;
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = ((maxXProbe - minXProbe) / 2) + minXProbe;
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = ((elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2) - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = ((maxZProbe - minZProbe) / 2) + minZProbe;
				paths[iLaw].z[1] = fbhValues[6][iLaw];
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
				if (maxMinElems - minElems[iElem] == NAN)
				{
					return PLUGIN_INVALID_ANGLE;
				}
                laws[iLaw].delays[iElem] = maxMinElems - minElems[iElem];
            }

			// Get the value of remarkable x,y,z coordinates when the number of element is odd.
			if (numberOfElements % 2 != 0) {
				paths[iLaw].x[0] = ((maxXProbe - minXProbe) / 2) + minXProbe;
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2] - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = ((maxZProbe - minZProbe) / 2) + minZProbe;
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

			// Get the value of remarkable x,y,z coordinates when the number of element is peer.
			else
			{
				paths[iLaw].x[0] = ((maxXProbe - minXProbe) / 2) + minXProbe;
				paths[iLaw].x[1] = fbhValues[5][iLaw];
				paths[iLaw].x[2] = xDef[iLaw];

				paths[iLaw].y[0] = ((elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2) - coupling.height;
				paths[iLaw].y[1] = fbhValues[4][iLaw];
				paths[iLaw].y[2] = yDef[iLaw];

				paths[iLaw].z[0] = ((maxZProbe - minZProbe) / 2) + minZProbe;
				paths[iLaw].z[1] = fbhValues[6][iLaw];
				paths[iLaw].z[2] = zDef[iLaw];
			}

            free(minElems);
        }
		#endif

		free(xDef);
		free(yDef);
		free(zDef);
		
		free(zIntB);
		
        for (int i = 0; i < fbhValues.size(); i++)
		{
			free(fbhValues[i]);
		}
	}

	// Retreive original coordinates
	for (int i = 0; i < numberOfElements; i++)
	{
		elements.coordinates.y[i] = -1.0 * elements.coordinates.y[i] + coupling.height;
	}
	
	

	return PLUGIN_NO_ERROR;
}

