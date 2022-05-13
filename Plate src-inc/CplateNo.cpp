#include <iostream>
#include <cmath>
#include <cstring>

#define PI 3.14159265358979323846

using namespace std;


/**
* maxArray
*
* get the max element from a vector
*
* @param vector<double> vector : The vector you want to get the max element
*
* @return double max : The max element of the vector
**/
double maxArray(const double* vector, int size) {
	double max = vector[0];

	for (int i = 0; i < size; i++) {
		if (max <= vector[i]) {
			max = vector[i];
		}
	}

	return max;
}


/**
* minArray
*
* get the min element from a vector
*
* @param vector<double> vector : The vector you want to get the min element
*
* @return double min : The min element of the vector
**/
double minArray(const double* vector, int size) {
	double min = vector[0];

	for (int i = 0; i < size; i++) {
		if (min > vector[i]) {
			min = vector[i];
		}
	}
	return min;
}


double* append(double* ar1, double* ar2, int len1, int len2) {
	double* arTmp = (double*)malloc(sizeof(double) * (len1 + len2));

	for (int i = 0; i < len1; i++) {
		arTmp[i] = ar1[i];
	}

	int cpt = 0;

	for (int i = len1; i < len1 + len2; i++) {
		arTmp[i] = ar2[cpt];
		cpt++;
	}

	return arTmp;
}



/**
* PlateDelayLawsCalculation
*
* compute all the delay laws for a plate
*
* @param int nbCoord : Length of xProbes array
* @param int nbLaw : Length of plates array
* @param double(*plates)[3] plates : Array of angle,obliquity,focus for all the angles
* @param int focalPoint : Focal point (US Path, Depth)
* @param double plateCelerity : Celerity of ultrasound in the plate
* @param double intCelerity : Celerity of ultrasound in the interface
* @param double intPath : Distance between the plate and the probe
* @param double* xProbes : Array of all probe elements X coordinate
* @param double* yProbes : Array of all probe elements Y coordinate
* @param double* zProbes : Array of all probe elements Z coordinate
*
* @return double** plateDelayLawss : Array of all delay laws for each probe element
**/
double** PlateDelayLawsCalculation(int nbCoord,
	int nbLaw,
	double* angle,
	double* obliquity,
	double* focus,
	int focalPoint,
	double plateCelerity,
	double intCelerity,
	double intPath,
	double* xProbe,
	double* yProbe,
	double* zProbe
)
{
	double** plateDelayLaws = (double**)malloc(sizeof(double) * nbLaw * nbCoord);


	// This condition handle the case when the distance between the probe and the plate equals to 0
	if (intPath == 0)
	{

		// For loop that convert obliquity, angle in radiant and if focalPoint equals to 1 we switch from US Path to Depth
		for (int iLaw = 0; iLaw < nbLaw; iLaw++) {
			obliquity[iLaw] = (obliquity[iLaw] / 180) * PI;
			angle[iLaw] = (angle[iLaw] / 180) * PI;


			switch (focalPoint) {
			case 1:
				focus[iLaw] = focus[iLaw] / cos(angle[iLaw]);
				break;

			case 0:
				break;
			}
		}


		// For loop that iter the number of delay law that we want
		double* distances;
		double* probeDelayLaws;

		// For loop that iter the number of delay law that we want
		for (int iLaw = 0; iLaw < nbLaw; iLaw++) {
			distances = (double*)malloc(nbCoord * sizeof(double));

			// For loop that compute all distances between probe elements and the focal point
			for (int j = 0; j < nbCoord; j++) {
				distances[j] = sqrt(pow((cos(angle[iLaw]) * focus[iLaw]) - xProbe[j], 2.0)
					+ pow((cos(obliquity[iLaw]) * (sin(angle[iLaw]) * focus[iLaw])) - yProbe[j], 2.0)
					+ pow(((obliquity[iLaw] * sin(angle[iLaw])) * sin(obliquity[iLaw])) - zProbe[j], 2.0));
			}

			
			plateDelayLaws[iLaw] = new double[nbCoord];
			double* delay = plateDelayLaws[iLaw];

			// For loop that push the delay law for each element of the probe
			for (int k = 0; k < nbCoord; k++) {
				delay[k] = ((maxArray(distances, nbCoord) - distances[k]) / plateCelerity) * 1000;
			}

			// TODO commenter
			free(distances);
		}
	}



	// This else handle the case when intPath is != 0
	else {

		// For loop that convert obliquity, angle in radiant and put them in arrays
		for (int p = 0; p < nbLaw; p++) {
			angle[p] = (angle[p] / 180) * PI;

			obliquity[p] = (obliquity[p] / 180) * PI;
		}

		double* xFbhP = (double*)malloc(nbLaw * sizeof(double));
		double* yFbhP = (double*)malloc(nbLaw * sizeof(double));
		double* zFbhP = (double*)malloc(nbLaw * sizeof(double));

		// For loop that compute all flat bottom hole and put them in arrays
		for (int i = 0; i < nbLaw; i++) {
			double xFhb = (cos(angle[i]) * focus[i]) + intPath;

			double yFhb = ((tan(asin((sin(angle[i]) * intCelerity) / plateCelerity)) * intPath)
				+ (sin(angle[i]) * focus[i]))
				* cos(obliquity[i]);

			double zFhb = ((tan(asin((sin(angle[i]) * intCelerity) / plateCelerity)) * intPath)
				+ (sin(angle[i]) * focus[i]))
				* sin(obliquity[i]);

			xFbhP[i] = xFhb;
			yFbhP[i] = yFhb;
			zFbhP[i] = zFhb;
		}

		// array1, array2, array3, array4 are where we are going to put all min / max boundaries for interest zones on the interface
		double* array1 = (double*)malloc(nbLaw * sizeof(double));
		double* array2 = (double*)malloc(nbLaw * sizeof(double));
		double* array3 = (double*)malloc(nbLaw * sizeof(double));
		double* array4 = (double*)malloc(nbLaw * sizeof(double));

		int cptAr1 = 0;
		int cptAr2 = 0;
		int cptAr3 = 0;
		int cptAr4 = 0;


		// For loop that determines min / max boundaries for the interest zone on the interface
		for (int i = 0; i < nbLaw; i++) {
			if (maxArray(yProbe, nbCoord) < yFbhP[i]) {
				array1[i] = ((intPath - minArray(xProbe, nbCoord))
					* ((yFbhP[i] - maxArray(yProbe, nbCoord))
						/ (xFbhP[i] - minArray(xProbe, nbCoord))))
					+ maxArray(yProbe, nbCoord);
				cptAr1++;
			}
			else {
				array1[i] = maxArray(yProbe, nbCoord);
				cptAr1++;
			}

			if (yFbhP[i] < minArray(yProbe, nbCoord)) {
				array2[i] = ((minArray(yProbe, nbCoord) + yFbhP[i])
					/ (xFbhP[i] - minArray(xProbe, nbCoord)))
					* (intPath - minArray(xProbe, nbCoord))
					+ minArray(yProbe, nbCoord);
				cptAr2++;
			}
			else {
				array2[i] = minArray(yProbe, nbCoord);
				cptAr2++;
			}

			if (zFbhP[i] > maxArray(zProbe, nbCoord)) {
				array3[i] = (maxArray(zProbe, nbCoord)
					+ ((intPath - minArray(xProbe, nbCoord))
						* ((zFbhP[i] - maxArray(zProbe, nbCoord))
							/ (xFbhP[i] - minArray(xProbe, nbCoord)))));
				cptAr3++;
			}
			else {
				array3[i] = maxArray(zProbe, nbCoord);
				cptAr3++;
			}

			if (zFbhP[i] < minArray(zProbe, nbCoord)) {
				array4[i] = (((zFbhP[i] + minArray(zProbe, nbCoord))
					/ (xFbhP[i] - minArray(xProbe, nbCoord)))
					* (intPath - minArray(xProbe, nbCoord)))
					+ minArray(zProbe, nbCoord);
				cptAr4++;
			}
			else {
				array4[i] = minArray(zProbe, nbCoord);
				cptAr4++;
			}
		}


		double subArray1 = maxArray(array1, cptAr1) - minArray(array2, cptAr2); // Size of the interest zone on the main axis
		double subArray2 = maxArray(array3, cptAr3) - minArray(array4, cptAr4); // Size of the interest zone on the another axis

		cout << subArray1 << endl;
		cout << subArray2 << endl;

		free(array1);
		free(array3);

		double resolution = 0.1;


		double* arrayIntPath = (double*)malloc(((subArray1 / resolution) + 1 + 1) * sizeof(double)); // Array where we are going to put the interface path
		double* arrayPreYInterestP = (double*)malloc(((subArray1 / resolution) + 1 + 1) * sizeof(double)); // Array that contain Y coodinates of interest points


		// For loop that compute X,Y interface point coordinate
		for (int i = 0; i < ((subArray1 / resolution) + 1); i++) {
			arrayIntPath[i] = intPath;

			arrayPreYInterestP[i] = (resolution * i) + minArray(array2, cptAr2);
		}

		free(array2);

		double* xIntP = (double*)malloc((((subArray1 / resolution) + 1 + 1) * ((subArray2 / resolution) + 1)) * sizeof(double));
		double* yIntP = (double*)malloc((((subArray1 / resolution) + 1 + 1) * ((subArray2 / resolution) + 1)) * sizeof(double));
		double* zIntP = (double*)malloc((((subArray1 / resolution) + 1 + 1) * ((subArray2 / resolution) + 1)) * sizeof(double));

		double* zIntPTmp = (double*)malloc(((subArray1 / resolution) + 1) * sizeof(double));

		// For loop that duplicate X,Y coordinates and compute Z interface point coordinate for matrix probe
		for (int i = 0; i < ((subArray2 / resolution) + 1); i++) {

			xIntP = append(xIntP, arrayIntPath, (int)(subArray1 / resolution) * (i), (subArray1 / resolution));
			yIntP = append(yIntP, arrayPreYInterestP, (int)(subArray1 / resolution) * (i), (subArray1 / resolution));


			for (int j = 0; j < ((subArray1 / resolution) + 1); j++) {
				zIntPTmp[j] = (i * resolution) + minArray(array4, cptAr4);
			}

			zIntP = append(zIntP, zIntPTmp, (int)(subArray1 / resolution) * (i), (subArray1 / resolution));
		}

		xIntP = append(xIntP, arrayIntPath, (int)(subArray1 / resolution) * ((subArray2 / resolution) + 1), (subArray1 / resolution));
		yIntP = append(yIntP, arrayPreYInterestP, (int)(subArray1 / resolution) * ((subArray2 / resolution) + 1), (subArray1 / resolution));
		zIntP = append(zIntP, zIntPTmp, (int)(subArray1 / resolution) * ((subArray2 / resolution) + 1), (subArray1 / resolution));

		(int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1);

		free(array4);
		free(arrayIntPath);
		free(arrayPreYInterestP);

		for (int iLaw = 0; iLaw < nbLaw; iLaw++) {

			double* timeFbhInt = (double*)malloc((int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the interest point to the flat bottom hole

			// For loop that compute the time for the US to go from the interest point to the flat bottom hole
			for (int j = 0; j < (int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1); j++) {
				timeFbhInt[j] = (sqrt(pow(xFbhP[iLaw] - xIntP[j], 2.0)
					+ pow(yFbhP[iLaw] - yIntP[j], 2.0)
					+ pow(zFbhP[iLaw] - zIntP[j], 2.0)))
					/ (plateCelerity / 1000);
				
			}


			double* delayLaw = (double*)malloc(nbCoord * sizeof(double)); // Array that contain delay law for each probe
			double* addTimeVector = (double*)malloc((int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the probe element to the flat bottom hole

			// For loop that compute the time for the US to go from the probe element to the interest point
			for (int probeElem = 0; probeElem < nbCoord; probeElem++) {

				double* timeProbeInt = (double*)malloc((int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1) * sizeof(double)); // Array that contain the time for the US to go from the probe element to the interest point

				// For loop that compute the time for the US to go from the probe element to the interest point
				for (int j = 0; j < (int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1); j++) {
					timeProbeInt[j] = (sqrt(pow(xProbe[probeElem] - xIntP[j], 2.0)
						+ pow(yProbe[probeElem] - yIntP[j], 2.0)
						+ pow(zProbe[probeElem] - zIntP[j], 2.0)))
						/ (intCelerity / 1000);
				}



				// For loop that add timeFbhInt and timeProbeInt
				for (int l = 0; l < (int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1); l++) {
					addTimeVector[l] = timeFbhInt[l] + timeProbeInt[l];
				}

                

				free(timeProbeInt);
               
				delayLaw[probeElem] = minArray(addTimeVector, (int)((subArray2 / resolution) + 1) * (int)((subArray1 / resolution) + 1));
			}

			free(timeFbhInt);
			free(addTimeVector);


			plateDelayLaws[iLaw] = new double[nbCoord];
			double* delay = plateDelayLaws[iLaw];

			// For loop that push the delay law for each element of the probe
			for (int dL = 0; dL < nbCoord; dL++) {
				delay[dL] = maxArray(delayLaw, nbCoord) - delayLaw[dL];
			}

			free(delayLaw);
		}

		free(xIntP);
		free(yIntP);
		free(zIntP);
		//free(zIntPTmp);

		free(xFbhP);
		free(yFbhP);
		free(zFbhP);

	}

	return plateDelayLaws;
}


int main() {
	double pla[16][3] = {{0,0,50}, {2,0,50}, {4,0,50}, {6,0,50}, {8,0,50}, {10,0,50}, {12,0,50}, {14,0,50}, {16,0,50}, {18,0,50},
	{20,0,50}, {22,0,50}, {24,0,50}, {26,0,50}, {28,0,50}, {30,0,50} };


	double ang[16] = { 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 };
	double obl[16] = { 0,0,0,0 ,0,0,0,0, 0,0,0,0, 0,0,0,0 };
	double foc[16] = { 50,50,50,50,50,50,50,50, 50,50,50,50, 50,50,50,50 };

	double yP[64]{
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
	0.3,
	0.9,
	1.5,
	2.1,
	2.7,
	3.3,
	3.9,
	4.5,
	5.1,
	5.7,
	6.3,
	6.9,
	7.5,
	8.1,
	8.7,
	9.3,
	9.9,
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
	18.9 };
	double xP[64] = { 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0 };
	double zP[64] = { 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0 };
	
	/*
	double ang[8] = { 0,2,4,6,8,10,12,14 };
	double obl[8] = { 0,0,0,0,0,0,0,0 };
	double foc[8] = { 50,50,50,50,50,50,50,50 };

	double xP[32] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	double yP[32] = { 3.15, 3.15, 3.15, 3.15, 2.25, 2.25, 2.25, 2.25, 1.35, 1.35, 1.35, 1.35, 0.45, 0.45, 0.45, 0.45, -0.45, -0.45, -0.45, -0.45, -1.35, -1.35,
					 -1.35, -1.35, -2.25, -2.25, -2.25, -2.25, -3.15, -3.15, -3.15, -3.15 };

	double zP[32] = { 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45,
				-0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45,-1.35 };
*/

	// Sonde sectorielle 128 éléments
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

*/

	clock_t timeReq;
	timeReq = clock();
	double** res = PlateDelayLawsCalculation(64, 16, ang, obl, foc, 0, 5940, 1500, 10, xP, yP, zP);
	timeReq = clock() - timeReq;
	/*
	for (int i = 0; i < 16; i++) {
		double* vect = res[i];
		std::cout << "---------------------" << endl;
		for (int j = 0; j < 64; j++) {
			cout << vect[j] << endl;
		}
		free(vect);
	}
	*/
	free(res);
	cout << ((float)timeReq / CLOCKS_PER_SEC) * 1000 << " ms" << endl;
}