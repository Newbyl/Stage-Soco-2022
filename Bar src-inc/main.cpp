#include <iostream>
#include "Cbar.h"
#include "unit.h"

using namespace std;

int main()
{
	// 32 éléments matricielle
	/*
	double tilts[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double notch[16] = {0, 2, 4, 6};
	double pos[16] = {50, 50, 50, 50};

	double xP[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	double yP[32] = {3.15, 3.15, 3.15, 3.15, 2.25, 2.25, 2.25, 2.25, 1.35, 1.35, 1.35, 1.35, 0.45, 0.45, 0.45, 0.45, -0.45, -0.45, -0.45, -0.45, -1.35, -1.35,
					-1.35, -1.35, -2.25, -2.25, -2.25, -2.25, -3.15, -3.15, -3.15, -3.15};

	double zP[32] = {1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45,
				-0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45,-1.35};

	int i[32] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4};
	*/

	// 64 éléments
	/*
	double tilts[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double notch[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double pos[16] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};

	double yP[64]{-18.9, -18.3, -17.7, -17.1, -16.5, -15.9, -15.3, -14.7, -14.1, -13.5, -12.9, -12.3, -11.7, -11.1, -10.5, -9.9, -9.3, -8.7, -8.1, -7.5, -6.9, -6.3, -5.7, -5.1, -4.5, -3.9, -3.3, -2.7,
				-2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1, 2.7, 3.3, 3.9, 4.5, 5.1, 5.7, 6.3, 6.9, 7.5, 8.1, 8.7, 9.3, 9.9, 10.5, 11.1, 11.7, 12.3, 12.9, 13.5, 14.1, 14.7, 15.3, 15.9, 16.5, 17.1, 17.7, 18.3, 18.9};
	double xP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double zP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	int i[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	*/

	// Sonde sectorielle 127 éléments
	
	double tilts[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double notch[16] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
	double pos[16] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};

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



	int nbElem = 127;
	int nbLaw = 16;
	double matcel = 3230.0;
	double ch = 70;
	int nbPtRemar = 3;
	int aglTy = 1;
	int defTy = 0;
	int proTy = 2;
	double barDiam = 404;

	double *la = (double *)malloc(nbLaw * nbElem * sizeof(double));
	double *xRemarq = (double *)malloc(3 * nbLaw * sizeof(double));
	double *yRemarq = (double *)malloc(3 * nbLaw * sizeof(double));
	double *zRemarq = (double *)malloc(3 * nbLaw * sizeof(double));

	Cbar cbar;
	cbar.Open();

	cbar.Set("Elements.coordinates.x", UNIT_mm, &nbElem, xP);
	cbar.Set("Elements.coordinates.y", UNIT_mm, &nbElem, yP);
	cbar.Set("Elements.Coordinates.z", UNIT_mm, &nbElem, zP);
	cbar.Set("Elements.coordinates.i", UNIT_mm, &nbElem, i);
	
	cbar.Set("Targets.Tilts", UNIT_deg, &nbLaw, tilts);
	cbar.Set("Targets.NotchesAngles", UNIT_deg, &nbLaw, notch);
	cbar.Set("Targets.Positions", UNIT_mm, &nbLaw, pos);
	cbar.Set("Material.Velocity", UNIT_mps, &matcel);

	cbar.Set("Defect_Type", UNIT_mm, &defTy);
	cbar.Set("Angle_Type", UNIT_mm, &aglTy);
	cbar.Set("Probe_Type", UNIT_mm, &proTy);

	cbar.Set("Coupling.Height", UNIT_mm, &ch);
	cbar.Set("BarDiameter", UNIT_mm, &barDiam);

	clock_t timeReq;
	timeReq = clock();
	cbar.ExecSync("Calculate");
	timeReq = clock() - timeReq;

	cbar.Get("laws", UNIT_ns, &nbLaw, &nbElem, la);

	cbar.Get("Paths.x", UNIT_mm, &nbLaw, &nbPtRemar, xRemarq);
	cbar.Get("Paths.y", UNIT_mm, &nbLaw, &nbPtRemar, yRemarq);
	cbar.Get("Paths.z", UNIT_mm, &nbLaw, &nbPtRemar, zRemarq);
	
	
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

	cbar.Close();

	free(yRemarq);
	free(xRemarq);
	free(zRemarq);
	free(la);
}
