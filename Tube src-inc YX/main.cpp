#include <iostream>
#include "Ctube.h"
#include "unit.h"

using namespace std;

int main()
{
	double yP[64]{-18.9, -18.3, -17.7, -17.1, -16.5, -15.9, -15.3, -14.7, -14.1, -13.5, -12.9, -12.3, -11.7, -11.1, -10.5, -9.9, -9.3, -8.7, -8.1, -7.5, -6.9, -6.3, -5.7, -5.1, -4.5, -3.9, -3.3, -2.7,
				  -2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1, 2.7, 3.3, 3.9, 4.5, 5.1, 5.7, 6.3, 6.9, 7.5, 8.1, 8.7, 9.3, 9.9, 10.5, 11.1, 11.7, 12.3, 12.9, 13.5, 14.1, 14.7, 15.3, 15.9, 16.5, 17.1, 17.7, 18.3, 18.9};

	double xP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double zP[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	double ang[16] = {18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18};
	double obl[16] = {0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5};

	int nbElem = 64;
	int nbLaw = 16;
	double ch = 35;
	int nbPtRemar = 3;
	double diam = 100;
	double thickness = 16;
	double focalL = 10000;


	double *la = (double *)malloc(nbLaw * nbElem * sizeof(double));

	double *alS = (double *)malloc(nbLaw * sizeof(double));
	double *xi = (double *)malloc(nbLaw * sizeof(double));
	double *yi = (double *)malloc(nbLaw * sizeof(double));
	double *zi = (double *)malloc(nbLaw * sizeof(double));





	Ctube ctube;

	ctube.Open();

	ctube.Set("Elements.coordinates.x", UNIT_mm, &nbElem, xP);

	ctube.Set("Elements.coordinates.y", UNIT_mm, &nbElem, yP);

	ctube.Set("Elements.coordinates.z", UNIT_mm, &nbElem, zP);

	ctube.Set("Targets.Tilts", UNIT_deg, &nbLaw, ang);
	ctube.Set("Targets.Skews", UNIT_deg, &nbLaw, obl);
	
	ctube.Set("diameter", UNIT_mm, &diam);
	ctube.Set("thickness", UNIT_mm, &thickness);
	ctube.Set("focal.length.coupling", UNIT_mm, &focalL);
	ctube.Set("coupling.height", UNIT_mm, &ch);

	clock_t timeReq;
	timeReq = clock();
	ctube.ExecSync("Calculate");
	timeReq = clock() - timeReq;

	ctube.Get("laws", UNIT_ns, &nbLaw, &nbElem, la);

	ctube.Get("alphaS", UNIT_deg, &nbLaw, alS);
	ctube.Get("xi3D", UNIT_mm, &nbLaw, xi);
	ctube.Get("yi3D", UNIT_mm, &nbLaw, yi);
	ctube.Get("zi3D", UNIT_mm, &nbLaw, zi);

	
	for (int i = 0; i < nbLaw * nbElem; i++)
	{
		cout << la[i] << endl;
	}

	/*
	for (int i = 0; i < nbLaw; i++)
	{
		cout << "AlphaS de la " << " loi " << i+1 << " : "<< alS[i] << endl;
	}
	*/
	
	for (int i = 0; i < nbLaw; i++)
	{
		cout << "xi3D de la " << " loi " << i+1 << " : "<< xi[i] << endl;
		cout << "yi3D de la " << " loi " << i+1 << " : "<< yi[i] << endl;
		cout << "zi3D de la " << " loi " << i+1 << " : " << zi[i] << endl;
	}
	

	cout << "Temps total ExecSync() + Calculate() : " 
	<< ((float)timeReq / CLOCKS_PER_SEC) * 1000 << " ms" << endl;
	

	free(alS);
	free(la);
	free(xi);
	free(yi);
	free(zi);

	ctube.Close();
}
