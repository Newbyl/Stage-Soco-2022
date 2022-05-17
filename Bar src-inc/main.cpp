#include <iostream>
#include "Cbar.h"
#include "unit.h"

using namespace std;

int main()
{
	// 32 éléments matricielle
	
	double tilts[2] = {0,2};
	double notch[2] = {0,2};
	double pos[2] = {1, 1};

	double xP[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	double yP[32] = {3.15, 3.15, 3.15, 3.15, 2.25, 2.25, 2.25, 2.25, 1.35, 1.35, 1.35, 1.35, 0.45, 0.45, 0.45, 0.45, -0.45, -0.45, -0.45, -0.45, -1.35, -1.35,
					-1.35, -1.35, -2.25, -2.25, -2.25, -2.25, -3.15, -3.15, -3.15, -3.15};

	double zP[32] = {1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45,
				-0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45,-1.35};

	int i[32] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4};


	int nbElem = 32;
	int nbLaw = 2;
	double matcel = 3230.0;
	double ch = 70;
	int nbPtRemar = 3;
	int aglTy = 0;
	int defTy = 0;
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
