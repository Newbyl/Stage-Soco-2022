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
double targetsNotcheAngle[2] = {0, 2};

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

        cout << z[iLaw] << endl;
	}

	std::vector<double*> values{x, y, z, utAngle};

	return values;
}

int main()
{
    fbhBuilder(barDiameter/2);
}