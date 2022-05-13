#include "framework.h"
#include "Ctube.h"
#include "unit.h"
#include "error_codes.h"
#include "framework.h"

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
 * @brief This fonction return all remarkable point of the US in the tube and the incidence angle.
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
	double xd3D = abs(m1PhC(0));
	double yd3D = abs(m1PhC(1));
	double zd3D = abs(m1PhC(2));
	
	// Coordinates of the reflex in the tube
	double xAlpha3D = abs(m2PhC(0));
	double yAlpha3D = abs(m2PhC(1));
	double zAlpha3D = abs(m2PhC(2));

	// Coordonné du point d'interface
	double xi3D = abs(m3PhC(0));
	double yi3D = abs(m3PhC(1));
	double zi3D = abs(m3PhC(2));

	// Angle of incidence on the surface.
	double alphaS = acos(xi3D / sqrt(pow(xi3D, 2.0) + pow(yi3D, 2.0) + pow(zi3D, 2.0))) / M_PI * 180;

	// Here we get all the values in a vector to use them in Calculate function.
	vector<double> returnVal {xi3D, yi3D, zi3D, alphaS, xd3D, yd3D, zd3D, xAlpha3D, yAlpha3D, zAlpha3D};

	return returnVal;
}


int Ctube::Calculate()
{
	using namespace std;
		
	// For loop that iterate number of target times.
	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
		std::vector<double> elipse = Ctube::newElipse(targets.skews[iLaw], targets.tilts[iLaw]);
		
		// Here we get all the values given by the function newElipse
		double xI3Dv = elipse.at(0);
		double yI3Dv = elipse.at(1);
		double zI3Dv = elipse.at(2);
		double alphaSv = elipse.at(3);

		double xD3Dv = elipse.at(4);
		double yD3Dv = elipse.at(5);
		double zD3Dv = elipse.at(6);

		double xR3Dv = elipse.at(7);
		double yR3Dv = elipse.at(8);
		double zR3Dv = elipse.at(9);


		double if1;
		double if2;

		if (targets.skews[iLaw] <= 360 && targets.skews[iLaw] >= 180)
			if2 = -(zI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
		else
			if2 = (zI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
		
		if (targets.skews[iLaw] <= 270 && targets.skews[iLaw] >= 90)
			if1 = -(yI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
		else
			if1 = (yI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0))));
		

		// Array of all the distances between the probe and the defect.
		double* distancesArray = (double*)malloc(numberOfElements * sizeof(double));

		
		// For loop where are going to compute all possible path.
		for (int iElem = 0; iElem < numberOfElements; iElem++)
		{
			double dist2 = sqrt(pow(((xI3Dv / (sqrt(pow(xI3Dv, 2.0) + pow(yI3Dv, 2.0) + pow(zI3Dv, 2.0)))) * focal.length.coupling)
			 - elements.coordinates.x[iElem], 2.0) 
			+ pow((if1 * focal.length.coupling) - elements.coordinates.y[iElem], 2.0)
			+ pow((if2 * focal.length.coupling) - elements.coordinates.z[iElem], 2.0));

			distancesArray[iElem] = dist2;
		}
		
		// Maximum value of tabDist array
		double maxDistancesArray = maxArray(distancesArray, numberOfElements);

		// For loop that iterate number of element times and where we are going to compute
		// all the delay for each element.
		for (int iElem = 0; iElem < numberOfElements; iElem++)
		{
			distancesArray[iElem] = (((maxDistancesArray - distancesArray[iElem]) / coupling.velocity) * 1000);

			laws[iLaw].delays[iElem] = distancesArray[iElem];
		}
		
		// Free of the memory used by distancesArray because we dont need it anymore.
		free(distancesArray);

		// Here we get the values from the newElipse function for each law.
		alphaS[iLaw] = alphaSv;
		xi3D[iLaw] = xI3Dv;
		yi3D[iLaw] = yI3Dv;
		zi3D[iLaw] = zI3Dv;
		
		// Get the value of remarkable x,y,z coordinates when the number of element is odd.
		if (numberOfElements % 2 != 0) {
			paths[iLaw].x[0] = elements.coordinates.x[numberOfElements / 2];
			paths[iLaw].x[1] = xI3Dv;
			paths[iLaw].x[2] = xR3Dv;
			paths[iLaw].x[3] = xD3Dv;

			paths[iLaw].y[0] = elements.coordinates.y[numberOfElements / 2];
			paths[iLaw].y[1] = yI3Dv;
			paths[iLaw].y[2] = yR3Dv;
			paths[iLaw].y[3] = yD3Dv;

			paths[iLaw].z[0] = elements.coordinates.z[numberOfElements / 2];
			paths[iLaw].z[1] = zI3Dv;
			paths[iLaw].z[2] = zR3Dv;
			paths[iLaw].z[3] = zD3Dv;
		}

		// Get the value of remarkable x,y,z coordinates when the number of element is peer.
		else
		{
			paths[iLaw].x[0] = (elements.coordinates.x[numberOfElements / 2] + elements.coordinates.x[(numberOfElements / 2) - 1]) / 2;
			paths[iLaw].x[1] = xI3Dv;
			paths[iLaw].x[2] = xR3Dv;
			paths[iLaw].x[3] = xD3Dv;

			paths[iLaw].y[0] = (elements.coordinates.y[numberOfElements / 2] + elements.coordinates.y[(numberOfElements / 2) - 1]) / 2;
			paths[iLaw].y[1] = yI3Dv;
			paths[iLaw].y[2] = yR3Dv;
			paths[iLaw].y[3] = yD3Dv;

			paths[iLaw].z[0] = (elements.coordinates.z[numberOfElements / 2] + elements.coordinates.z[(numberOfElements / 2) - 1]) / 2;
			paths[iLaw].z[1] = zI3Dv;
			paths[iLaw].z[2] = zR3Dv;
			paths[iLaw].z[3] = zD3Dv;
		}

	}

	return PLUGIN_NO_ERROR;
}

