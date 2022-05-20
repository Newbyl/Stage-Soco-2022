#include <iostream>
#include "Cplate.h"
#include "unit.h"

using namespace std;

int main()
{
	// 64 éléments
	/*
	double ang[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double obl[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double foc[16] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};

	double yP[64]{-18.9, -18.3, -17.7, -17.1, -16.5, -15.9, -15.3, -14.7, -14.1, -13.5, -12.9, -12.3, -11.7, -11.1, -10.5, -9.9, -9.3, -8.7, -8.1, -7.5, -6.9, -6.3, -5.7, -5.1, -4.5, -3.9, -3.3, -2.7,
				-2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1, 2.7, 3.3, 3.9, 4.5, 5.1, 5.7, 6.3, 6.9, 7.5, 8.1, 8.7, 9.3, 9.9, 10.5, 11.1, 11.7, 12.3, 12.9, 13.5, 14.1, 14.7, 15.3, 15.9, 16.5, 17.1, 17.7, 18.3, 18.9};

	double xP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double zP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int i[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	*/
	// 32 éléments matricielle
	
	double ang[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double obl[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double foc[16] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};

	double xP[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	double yP[32] = {3.15, 3.15, 3.15, 3.15, 2.25, 2.25, 2.25, 2.25, 1.35, 1.35, 1.35, 1.35, 0.45, 0.45, 0.45, 0.45, -0.45, -0.45, -0.45, -0.45, -1.35, -1.35,
					-1.35, -1.35, -2.25, -2.25, -2.25, -2.25, -3.15, -3.15, -3.15, -3.15};

	double zP[32] = {1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45,
				-0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45,-1.35};

	int i[32] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4};
	
	// 128 éléments
	/*
	double ang[32] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62};
	double obl[32] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double foc[32] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};

	double xP[128] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	double zP[128] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	int i[128] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

	double yP[128] = {
		-38.1,
	-37.5,
	-36.9,
	-36.3,
	-35.7,
	-35.1,
	-34.5,
	-33.9,
	-33.3,
	-32.7,
	-32.1,
	-31.5,
	-30.9,
	-30.3,
	-29.7,
	-29.1,
	-28.5,
	-27.9,
	-27.3,
	-26.7,
	-26.1,
	-25.5,
	-24.9,
	-24.3,
	-23.7,
	-23.1,
	-22.5,
	-21.9,
	-21.3,
	-20.7,
	-20.1,
	-19.5,
	-18.9,
	-18.3,
	-17.7,
	-17.1,
	-16.5,
	-15.9,
	-15.3,
	-14.7,
	-14.1,
	-13.5,
	-12.9,
	-12.3,
	-11.7,
	-11.1,
	-10.5,
	-9.9,
	-9.3,
	-8.7,
	-8.1,
	-7.5,
	-6.9,
	-6.3,
	-5.7,
	-5.1,
	-4.5,
	-3.9,
	-3.3,
	-2.7,
	-2.1,
	-1.5,
	-0.9,
	-0.3,
	0.3	,
	0.9	,
	1.5	,
	2.1	,
	2.7	,
	3.3	,
	3.9	,
	4.5	,
	5.1	,
	5.7	,
	6.3	,
	6.9	,
	7.5	,
	8.1	,
	8.7	,
	9.3	,
	9.9	,
	10.5,
	11.1,
	11.7,
	12.3,
	12.9,
	13.5,
	14.1,
	14.7,
	15.3,
	15.9,
	16.5,
	17.1,
	17.7,
	18.3,
	18.9,
	19.5,
	20.1,
	20.7,
	21.3,
	21.9,
	22.5,
	23.1,
	23.7,
	24.3,
	24.9,
	25.5,
	26.1,
	26.7,
	27.3,
	27.9,
	28.5,
	29.1,
	29.7,
	30.3,
	30.9,
	31.5,
	32.1,
	32.7,
	33.3,
	33.9,
	34.5,
	35.1,
	35.7,
	36.3,
	36.9,
	37.5,
	38.1
	};

*/
// Sonde sectorielle 127 éléments
/*
	double ang[8] = {0, 2, 4, 6, 8, 10, 12, 14};
	double obl[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	double foc[8] = {50, 50, 50, 50, 50, 50, 50, 50};

	double yP[127]{0,
0,
0.779423,
0.779423,
1.10218E-16,
-0.779423,
-0.779423,
0,
0.9,
1.55885,
1.8,
1.55885,
0.9,
2.20436E-16,
-0.9,
-1.55885,
-1.8,
-1.55885,
-0.9,
0,
0.923454,
1.73553,
2.33827,
2.65898,
2.65898,
2.33827,
1.73553,
0.923454,
3.30655E-16,
-0.923454,
-1.73553,
-2.33827,
-2.65898,
-2.65898,
-2.33827,
-1.73553,
-0.923454,
0,
0.931749,
1.8,
2.54558,
3.11769,
3.47733,
3.6,
3.47733,
3.11769,
2.54558,
1.8,
0.931749,
4.40873E-16,
-0.931749,
-1.8,
-2.54558,
-3.11769,
-3.47733,
-3.6,
-3.47733,
-3.11769,
-2.54558,
-1.8,
-0.931749,
0,
0.935603,
1.83031,
2.64503,
3.34415,
3.89711,
4.27975,
4.47535,
4.47535,
4.27975,
3.89711,
3.34415,
2.64503,
1.83031,
0.935603,
5.51091E-16,
-0.935603,
-1.83031,
-2.64503,
-3.34415,
-3.89711,
-4.27975,
-4.47535,
-4.47535,
-4.27975,
-3.89711,
-3.34415,
-2.64503,
-1.83031,
-0.935603,
0,
0.9377,
1.84691,
2.7,
3.47105,
4.13664,
4.67654,
5.07434,
5.31796,
5.4,
5.31796,
5.07434,
4.67654,
4.13664,
3.47105,
2.7,
1.84691,
0.9377,
6.61309E-16,
-0.9377,
-1.84691,
-2.7,
-3.47105,
-4.13664,
-4.67654,
-5.07434,
-5.31796,
-5.4,
-5.31796,
-5.07434,
-4.67654,
-4.13664,
-3.47105,
-2.7,
-1.84691,
-0.9377 
};

	double xP[127] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double zP[127] = {0,
0.9,
0.45,
-0.45,
-0.9,
-0.45,
0.45,
1.8,
1.55885,
0.9,
1.10218E-16,
-0.9,
-1.55885,
-1.8,
-1.55885,
-0.9,
-3.30655E-16,
0.9,
1.55885,
2.7,
2.53717,
2.06832,
1.35,
0.46885,
-0.46885,
-1.35,
-2.06832,
-2.53717,
-2.7,
-2.53717,
-2.06832,
-1.35,
-0.46885,
0.46885,
1.35,
2.06832,
2.53717,
3.6,
3.47733,
3.11769,
2.54558,
1.8,
0.931749,
2.20436E-16,
-0.931749,
-1.8,
-2.54558,
-3.11769,
-3.47733,
-3.6,
-3.47733,
-3.11769,
-2.54558,
-1.8,
-0.931749,
-6.61309E-16,
0.931749,
1.8,
2.54558,
3.11769,
3.47733,
4.5,
4.40166,
4.11095,
3.64058,
3.01109,
2.25,
1.39058,
0.470378,
-0.470378,
-1.39058,
-2.25,
-3.01109,
-3.64058,
-4.11095,
-4.40166,
-4.5,
-4.40166,
-4.11095,
-3.64058,
-3.01109,
-2.25,
-1.39058,
-0.470378,
0.470378,
1.39058,
2.25,
3.01109,
3.64058,
4.11095,
4.40166,
5.4,
5.31796,
5.07434,
4.67654,
4.13664,
3.47105,
2.7,
1.84691,
0.9377,
3.30655E-16,
-0.9377,
-1.84691,
-2.7,
-3.47105,
-4.13664,
-4.67654,
-5.07434,
-5.31796,
-5.4,
-5.31796,
-5.07434,
-4.67654,
-4.13664,
-3.47105,
-2.7,
-1.84691,
-0.9377,
-9.91964E-16,
0.9377,
1.84691,
2.7,
3.47105,
4.13664,
4.67654,
5.07434,
5.31796	
};

int i[127] = {1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
*/



	int nbElem = 32;
	int nbLaw = 16;
	double matcel = 5940.0;
	double ch = 10;
	int nbPtRemar = 3;
	int probT = 1;

	double *la = (double *)malloc(nbLaw * nbElem * sizeof(double));
	double *xRemarq = (double *)malloc(3 * nbLaw * sizeof(double));
	double *yRemarq = (double *)malloc(3 * nbLaw * sizeof(double));
	double *zRemarq = (double *)malloc(3 * nbLaw * sizeof(double));

	Cplate cplate;
	cplate.Open();

	cplate.Set("Elements.coordinates.x", UNIT_mm, &nbElem, yP);
	cplate.Set("Elements.coordinates.y", UNIT_mm, &nbElem, xP);
	cplate.Set("Elements.coordinates.z", UNIT_mm, &nbElem, zP);
	cplate.Set("Elements.coordinates.i", UNIT_mm, &nbElem, i);
	
	cplate.Set("Targets.Tilts", UNIT_deg, &nbLaw, ang);
	cplate.Set("Targets.Skews", UNIT_deg, &nbLaw, obl);
	cplate.Set("Targets.Positions", UNIT_mm, &nbLaw, foc);
	cplate.Set("Material.Velocity", UNIT_mps, &matcel);

	cplate.Set("probeType", UNIT_mm, &probT);
	cplate.Set("Coupling.Height", UNIT_mm, &ch);

	clock_t timeReq;
	timeReq = clock();
	cplate.ExecSync("Calculate");
	timeReq = clock() - timeReq;

	cplate.Get("laws", UNIT_ns, &nbLaw, &nbElem, la);

	cplate.Get("Paths.x", UNIT_mm, &nbLaw, &nbPtRemar, xRemarq);
	cplate.Get("Paths.y", UNIT_mm, &nbLaw, &nbPtRemar, yRemarq);
	cplate.Get("Paths.z", UNIT_mm, &nbLaw, &nbPtRemar, zRemarq);
	
	
	for (int i = 0; i < nbLaw * nbElem; i++)
	{
		cout << la[i] << endl;
	}
	
	/*
	for (int i = 2; i < nbLaw * 3; i+=3)
	{
		cout << "x : " << xRemarq[i] << endl;
		cout << "y : " << yRemarq[i] << endl;
	}
	*/
	/*
	for (int i = 0; i < nbLaw * 3 ; i++)
	{
		//cout << "x : " << xRemarq[i] << endl;
		cout << "x : " << xRemarq[i] << endl;
	}
	*/

	

	cout << "temps total : " << ((float)timeReq / CLOCKS_PER_SEC) * 1000 << " ms" << endl;

	cplate.Close();

	free(yRemarq);
	free(xRemarq);
	free(zRemarq);
	free(la);
}
