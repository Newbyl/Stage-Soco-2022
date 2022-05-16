#include <vector>
#include <iostream>
#include <cmath>

double numberOfTargets = 2;
int angleType = 1;
double couplingVelocity = 1480.0;
double materialVelocity = 3230;
double barDiameter = 404;
double couplingHeight = 70;
double targetsPositions[2] = {1, 1};
double targetsTilts[2] = {0, 2};
double targetsNotchesAngles[2] = {0, 2};

using namespace std;

std::vector<double*> fbhBuilder(double barDiameter2)
{
	double* ai = (double*)malloc(numberOfTargets * sizeof(double));
	double* ar = (double*)malloc(numberOfTargets * sizeof(double));

	double* utAngle = (double*)malloc(numberOfTargets * sizeof(double));

	double* x = (double*)malloc(numberOfTargets * sizeof(double));
	double* y = (double*)malloc(numberOfTargets * sizeof(double));
	double* z = (double*)malloc(numberOfTargets * sizeof(double));

	for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
	{
        ai[iLaw] = asin((couplingVelocity / materialVelocity) * sin(targetsTilts[iLaw] / 180 * M_PI));
        ar[iLaw] = targetsTilts[iLaw] / 180 * M_PI;
        
        

		double a = M_PI - ai[iLaw];
		double ad = a + asin(sin(a) * (barDiameter2 / (couplingHeight + barDiameter2)));

		double y1 = sin((M_PI - ad)) * barDiameter2;
		double x0 = cos((M_PI - ad)) * barDiameter2;
		double x1 = (barDiameter2 - x0);

		double x2 = sin((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barDiameter2);
		double y2 = cos((M_PI - (ar[iLaw] + ((M_PI / 2) - (M_PI - ad))))) * ((sin(M_PI - (ar[iLaw] * 2)) / sin(ar[iLaw])) * barDiameter2);

		double distance = sqrt(pow(x2, 2.0) + pow(y2, 2.0) + pow(0, 2.0));

		if (x0 == barDiameter2)
		{
			x[iLaw] = targetsPositions[iLaw];
			y[iLaw] = 0;
		}
		else
		{
			x[iLaw] = (targetsPositions[iLaw] * (x2 / distance)) + (barDiameter2 - x0);
			y[iLaw] = (targetsPositions[iLaw] * (y2 / distance)) + y1;
		}

		z[iLaw] = targetsPositions[iLaw] * (0 / distance);

		utAngle[iLaw] = ar[iLaw] / M_PI * 180;

        cout << utAngle[iLaw] << endl;
	}

	std::vector<double*> values{x, y, z, utAngle};

	return values;
}

std::vector<double*> notcheBuilder(double barDiameter2)
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
		if (angleType == 1)
		{
			ai[iLaw] = asin((couplingVelocity / materialVelocity) * sin(targetsNotchesAngles[iLaw] / 180 * M_PI));
			ar[iLaw] = targetsNotchesAngles[iLaw] / 180 * M_PI;
		}

		else if (angleType == 0)
		{
			ar[iLaw] = asin((materialVelocity / couplingVelocity) * sin(targetsNotchesAngles[iLaw] / 180 * M_PI));
			ai[iLaw] = targetsNotchesAngles[iLaw] / 180 * M_PI;
		}
		
		

		double a = M_PI - ai[iLaw];
		double ad = asin(sin(a) * (barDiameter2 / (couplingHeight + barDiameter2)));

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

		x[iLaw] = x2 + x1;
		y[iLaw] = y2 + y1;
		z[iLaw] = 0;
	}



	std::vector<double*> values{x, y, z, utAngle, focalLength};

	return values;
}

int main()
{
    // fbhBuilder(barDiameter/2);
    notcheBuilder(barDiameter/2);
}