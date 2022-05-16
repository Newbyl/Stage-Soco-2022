#include <vector>
#include <iostream>
#include <cmath>


int numberOfTargets = 2;
int numberOfElements = 32;
int angleType = 1;
double couplingVelocity = 1480.0;
double materialVelocity = 3230;
double barDiameter = 404;
double couplingHeight = 70;
double targetsPositions[2] = {1, 1};
double targetsTilts[2] = {0, 2};
double targetsNotchesAngles[2] = {0, 2};
int defectType = 0;
double resolution = 0.1;

// 32 éléments matricielle
	
double ang[8] = {0,2,4,6,8,10,12,14};
double obl[8] = {0,0,0,0,0,0,0,0};
double foc[8] = {50,50,50,50,50,50,50,50};

double xP[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double yP[32] = {3.15, 3.15, 3.15, 3.15, 2.25, 2.25, 2.25, 2.25, 1.35, 1.35, 1.35, 1.35, 0.45, 0.45, 0.45, 0.45, -0.45, -0.45, -0.45, -0.45, -1.35, -1.35,
                -1.35, -1.35, -2.25, -2.25, -2.25, -2.25, -3.15, -3.15, -3.15, -3.15};

double zP[32] = {1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45,
            -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45, -1.35, 1.35, 0.45, -0.45,-1.35};




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


double maxArray(double *array, int size)
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

double minArray(double *array, int size)
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


double *append(vector<double> ar1, double* ar2, int len1, int len2)
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

	//free(ar2);

	return arTmp;
}


void Calculate()
{
    double* delayLaws = (double*)malloc(numberOfTargets * numberOfElements * sizeof(double));

    if (defectType == 0){
        double* asinTiltRad = (double*)malloc(numberOfTargets * sizeof(double));
        double* zDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* xDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* yDef = (double*)malloc(numberOfTargets * sizeof(double));

        for (int i = 0; i < numberOfTargets; i++)
        {
            asinTiltRad[i] = asin(sin(targetsTilts[i] / 180 * M_PI) * (couplingVelocity / materialVelocity));
        }

        double maxAngle = maxArray(asinTiltRad, numberOfTargets);
        double minAngle = minArray(asinTiltRad, numberOfTargets);

        vector<double*> fbhValues = fbhBuilder(barDiameter/2);
        
        for (int i = 0; i < numberOfTargets; i++)
        {
            zDef[i] = fbhValues[2][i];
            xDef[i] = fbhValues[0][i] + couplingHeight;
            yDef[i] = fbhValues[1][i];
        }

        // A copier coller a partir d'ici pour le else.

        double if1;
        double if2;

        if (tan(maxAngle) * couplingHeight == tan(minAngle) * couplingHeight)
        {
            if (tan(maxAngle) * couplingHeight >= 0)
            {
                if1 = tan(maxAngle) * couplingHeight;
                if2 = 0;
            }
            else
            {
                if1 = 0;
                if2 = tan(minAngle) * couplingHeight;
            }
        }
        else
        {
            if1 = tan(maxAngle) * couplingHeight;
            if2 = tan(minAngle) * couplingHeight;
        }

        vector<double> preXIntB;
        vector<double> preYIntB;

        double minYProbe = minArray(yP, numberOfElements);
        double maxYProbe = maxArray(yP, numberOfElements);
        double minZProbe = minArray(zP, numberOfElements);
        double maxZProbe = maxArray(zP, numberOfElements);

        if (numberOfElements <= 1)
            if1 = resolution;
        if (numberOfElements > 1)
            if2 = abs(if2);
        else
            if2 = resolution;


        for (int i = 0; i < ((barDiameter * M_PI / 2) / resolution) + 1; i++)
        {
            if (cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                >= minYProbe - if2 &&
                cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                < maxYProbe + if1)
            {
                preXIntB.push_back(sin(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2) + (barDiameter / 2) + couplingHeight);

                preYIntB.push_back(cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2));

    
            }
        }

        double* xIntB;
        double* yIntB;
        double* zIntB = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
        * sizeof(double));



        for (int i = 0; i < (((maxZProbe - minZProbe) / resolution) + 1); i++)
        {
            xIntB = append(preXIntB, xIntB, preXIntB.size(), i * preXIntB.size());
            yIntB = append(preYIntB, yIntB, preYIntB.size(), i * preYIntB.size());

            for (int j = i * (preYIntB.size()); j < preYIntB.size() * (i + 1); j++)
            {
                zIntB[j] = (resolution * i) + minZProbe;
            }
        }

        for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
            * sizeof(double));

            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (materialVelocity / 1000);
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
                    distIntProbe[iIntPoint] = sqrt(pow(xP[iElem] - xIntB[iIntPoint], 2.0) 
                    + pow(yP[iElem] - yIntB[iIntPoint], 2.0) 
                    + pow(zP[iElem] - zIntB[iIntPoint], 2.0)) / (couplingVelocity / 1000);

                    addDist[iIntPoint] = distIntProbe[iIntPoint] + distDefInt[iIntPoint];
                }

                minElems[iElem] = minArray(addDist, (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size());
            }

            double maxMinElems = maxArray(minElems, numberOfElements);

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                delayLaws[iElem] = maxMinElems - minElems[iElem];
            }

            free(minElems);
        }
        
    }

    else
    {
        double* asinNotcheRad = (double*)malloc(numberOfTargets * sizeof(double));
        double* zDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* xDef = (double*)malloc(numberOfTargets * sizeof(double));
        double* yDef = (double*)malloc(numberOfTargets * sizeof(double));

        for (int i = 0; i < numberOfTargets; i++)
        {
            if (angleType == 1)
            {
                asinNotcheRad[i] = asin(sin(targetsNotchesAngles[i] / 180 * M_PI) * (couplingVelocity / materialVelocity));
            }
            else
            {
                asinNotcheRad[i] =targetsNotchesAngles[i] / 180 * M_PI;
            }
        }

        double maxAngle = maxArray(asinNotcheRad, numberOfTargets);
        double minAngle = minArray(asinNotcheRad, numberOfTargets);

        vector<double*> notcheValues = notcheBuilder(barDiameter/2);

        for (int i = 0; i < numberOfTargets; i++)
        {
            zDef[i] = notcheValues[2][i];
            xDef[i] = notcheValues[0][i] + couplingHeight;
            yDef[i] = notcheValues[1][i];
        }

        double* deflexionAngle = (double*)malloc(numberOfTargets * sizeof(double));

        for (int i = 0; i < numberOfTargets; i++)
        {
            if (angleType == 1)
            {
                deflexionAngle[i] = asin((sin(targetsNotchesAngles[i] / 180 * M_PI) / materialVelocity) * couplingVelocity);
            }
            else
            {
                deflexionAngle[i] = targetsNotchesAngles[i] / 180 * M_PI;
            }

            deflexionAngle[i] = round(asin(sin(M_PI - deflexionAngle[i]) * ((barDiameter / 2) / ((barDiameter / 2) + couplingHeight)))
                                / M_PI * 180);
        }

        

        double if1;
        double if2;

        if (tan(maxAngle) * couplingHeight == tan(minAngle) * couplingHeight)
        {
            if (tan(maxAngle) * couplingHeight >= 0)
            {
                if1 = tan(maxAngle) * couplingHeight;
                if2 = 0;
            }
            else
            {
                if1 = 0;
                if2 = tan(minAngle) * couplingHeight;
            }
        }
        else
        {
            if1 = tan(maxAngle) * couplingHeight;
            if2 = tan(minAngle) * couplingHeight;
        }

        vector<double> preXIntB;
        vector<double> preYIntB;

        double minYProbe = minArray(yP, numberOfElements);
        double maxYProbe = maxArray(yP, numberOfElements);
        double minZProbe = minArray(zP, numberOfElements);
        double maxZProbe = maxArray(zP, numberOfElements);

        if (numberOfElements <= 1)
            if1 = resolution;
        if (numberOfElements > 1)
            if2 = abs(if2);
        else
            if2 = resolution;


        for (int i = 0; i < ((barDiameter * M_PI / 2) / resolution) + 1; i++)
        {
            if (cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                >= minYProbe - if2 &&
                cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) * (barDiameter / 2)
                < maxYProbe + if1)
            {
                preXIntB.push_back(sin(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2) + (barDiameter / 2) + couplingHeight);

                preYIntB.push_back(cos(((M_PI / ((barDiameter * M_PI / 2) / resolution)) * i) - M_PI) 
                * (barDiameter / 2));

    
            }
        }

        double* xIntB;
        double* yIntB;
        double* zIntB = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
        * sizeof(double));



        for (int i = 0; i < (((maxZProbe - minZProbe) / resolution) + 1); i++)
        {
            xIntB = append(preXIntB, xIntB, preXIntB.size(), i * preXIntB.size());
            yIntB = append(preYIntB, yIntB, preYIntB.size(), i * preYIntB.size());

            for (int j = i * (preYIntB.size()); j < preYIntB.size() * (i + 1); j++)
            {
                zIntB[j] = (resolution * i) + minZProbe;
            }
        }

        for (int iLaw = 0; iLaw < numberOfTargets; iLaw++)
        {
            double* distDefInt = (double*)malloc((((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size() 
            * sizeof(double));

            for (int iIntPoint = 0; iIntPoint < (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size(); iIntPoint++)
            {
                distDefInt[iIntPoint] = sqrt(pow(xDef[iLaw] - xIntB[iIntPoint], 2.0) 
                + pow(yDef[iLaw] - yIntB[iIntPoint], 2.0) 
                + pow(zDef[iLaw] - zIntB[iIntPoint], 2.0)) / (materialVelocity / 1000);
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
                    distIntProbe[iIntPoint] = sqrt(pow(xP[iElem] - xIntB[iIntPoint], 2.0) 
                    + pow(yP[iElem] - yIntB[iIntPoint], 2.0) 
                    + pow(zP[iElem] - zIntB[iIntPoint], 2.0)) / (couplingVelocity / 1000);

                    addDist[iIntPoint] = distIntProbe[iIntPoint] + distDefInt[iIntPoint];
                }

                minElems[iElem] = minArray(addDist, (((maxZProbe - minZProbe) / resolution) + 1) * preYIntB.size());
                
            }

            double maxMinElems = maxArray(minElems, numberOfElements);

            for (int iElem = 0; iElem < numberOfElements; iElem++)
            {
                delayLaws[iElem] = maxMinElems - minElems[iElem];
            }

            free(minElems);
        }


    }
}


int main()
{
    // fbhBuilder(barDiameter/2);
    // notcheBuilder(barDiameter/2);

    clock_t timeReq;
	timeReq = clock();

    Calculate();

    timeReq = clock() - timeReq;
    cout << "temps total : " << ((float)timeReq / CLOCKS_PER_SEC) * 1000 << " ms" << endl;
}