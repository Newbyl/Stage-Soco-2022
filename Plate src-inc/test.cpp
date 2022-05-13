#include <cmath>
#include <iostream>
#include <vector>
#include "eigen-3.4.0/Eigen/LU"
#include "eigen-3.4.0/Eigen/Core"


double materialV = 3230;
double couplingV = 1500;
double couplingH = 35;
double diameter = 100;
double thickness = 16;

std::vector<double> newElipse(double skew, double alphaI)
{
	using namespace std;
	using namespace Eigen;

	// Faire attention aussi aux valeurs ecrit en dure ici a changer dans l'impl
	// Normalement ok.
	double penetrationAngle;
	double alphaR = 90.0 - asin(((materialV / couplingV) * sin(alphaI / 180.0 * M_PI))) / M_PI * 180.0;

	penetrationAngle = alphaR;
	

	double cga = cos(-M_PI/2.0);
	double sga = sin(-M_PI/2.0);	

	MatrixXf mat(4,4);
	mat << cga, -sga, 0, 0 ,
	sga, cga, 0, 0,
	0, 0 ,1 ,0,
	0, 0, 0, 1;


	double alphaObl = skew / 180.0 * M_PI;
	double betaProf = -((penetrationAngle / 180.0) * M_PI);

	MatrixXf dirD1D2(4,1);
	dirD1D2 << cos(alphaObl) * cos(betaProf), sin(betaProf), cos(betaProf) * sin(alphaObl), 0.0;

	MatrixXf dir2;
	dir2 = mat * dirD1D2;

	double teta = (atan2(dir2(1), dir2(0)) * 180.0) / M_PI;
	double alpha = (atan2(dir2(2), dir2(0)) * 180.0) / M_PI;

	
	double r1 = diameter / 2.0;
	double r2 = r1 - thickness;

	double xm1 = r1;
	double ym1 = 0.0;
	double zm1 = 0.0;

	MatrixXf m1(4,1);
	m1 << r1, 0.0, 0.0, 1.0;

	alpha = alpha / 180.0 * M_PI;
	teta = (teta / 180.0) * M_PI;

	// Si probleme faut revoir le sin(alpha) + bas parce que ca vient surement de la.
	double beta = atan(cos(alpha) * tan(teta));
	double b1 = sin(beta);
	double a1 = cos(alpha) * cos(beta);
	double c1 = cos(beta) * std::sin(alpha);

	double discri = sqrt(pow((xm1 * 2) * a1, 2.0) - ((pow(a1, 2.0) + pow(b1, 2.0)) * (pow(r1, 2.0) - pow(r2, 2.0)) * 4)) + -((xm1 * 2) * a1);

	double tPp = discri / ((pow(a1, 2.0) + pow(b1, 2.0)) * 2);
	double tP = (-((xm1 * 2.0) * a1) - 
	sqrt(pow((xm1 * 2) * a1, 2.0) - ((pow(a1, 2.0) + pow(b1, 2.0)) * 
	(pow(r1, 2.0) - pow(r2, 2.0)) * 4))) / 
	((pow(a1, 2.0) + pow(b1, 2.0)) * 2);

	double T_LgM1M2;

	if (tPp < tP)
		T_LgM1M2 = tPp;
	else
		T_LgM1M2 = tP;
	
	double xm2 = (T_LgM1M2 * a1) + xm1;
	double ym2 = b1 * T_LgM1M2;
	double zm2 = T_LgM1M2 * c1;

	MatrixXf m2(4,1);
	m2 << xm2, ym2, zm2, 1;


	double gamma2 = atan(ym2 / xm2);

	double sgam = sin(gamma2);
	double cgam = cos(gamma2);

	Matrix4f mat1;
	mat1 << cgam, sgam, 0, -r2,
	-sgam, cgam, 0, 0,
	0, 0, 1.0, -zm2,
	0, 0, 0, 1.0;

	MatrixXf m1p;
	m1p = mat1 * m1;

	
	MatrixXf m1pr(4, 1);
	m1pr << m1p(0), -m1p(1), -m1p(2), 1;

	MatrixXf m1r;
	m1r = mat1.inverse() * m1pr;

	double lg_m2_m1r = sqrt(pow(xm2 - m1r(0), 2.0) + pow(ym2 - m1r(1), 2.0) + pow(zm2 - m1r(2), 2.0));

	double a2 = (xm2 - m1r(0)) / lg_m2_m1r;
	double b2 = (ym2 - m1r(1)) / lg_m2_m1r;
	double c2 = (zm2 - m1r(2)) / lg_m2_m1r;

	double add1 = pow(xm2, 2.0) + pow(ym2, 2.0);

	double discriminant = (pow((ym2 * 2 * b2) + (xm2 * 2 * a2),2.0)) - 
	((add1 - pow(xm1, 2.0)) * (pow(a2, 2.0) + pow(b2, 2.0)) * 4);

	double tPp2 = (sqrt(discriminant) + -((ym2 * 2 * b2) + (xm2 * 2 * a2))) / ((pow(a2, 2.0) + pow(b2, 2.0)) * 2);

	double tP2 = (-((ym2 * 2 * b2) + (xm2 * 2 * a2)) - sqrt(discriminant)) / ((pow(a2, 2.0) + pow(b2, 2.0)) * 2);

	double t_lg_m2m3;

	if (tPp2 < tP2)
		t_lg_m2m3 = tPp2;
	else
		t_lg_m2m3 = tP2;

	double xm3 = a2 * t_lg_m2m3 + xm2;
	double ym3 = b2 * t_lg_m2m3 + ym2;
	double zm3 = c2 * t_lg_m2m3 + zm2;

	MatrixXf m3(4, 1);
	m3 << xm3, ym3, zm3, 1;

	double gamma3 = atan(ym3 / xm3);

	double sgam2 = sin(gamma3);
	double cgam2 = cos(gamma3);

	Matrix4f mat3;
	mat3 << cgam2, sgam2, 0, -xm1,
		-sgam2, cgam2, 0, 0,
		0, 0, 1, -zm3,
		0, 0, 0, 1;

	MatrixXf m2s;
	m2s = mat3 * m2;

	double lgM2M3 = sqrt(pow(m2s(0), 2.0) + pow(m2s(1), 2.0) + pow(m2s(2), 2.0));
	double lget = sqrt(pow(m2s(1), 2.0) + pow(m2s(2), 2.0));
	
	MatrixXf dirM2M3(4, 1);
	dirM2M3 << -m2s(0) / lgM2M3, -m2s(1) / lgM2M3, -m2s(2) / lgM2M3, 0;

	double aglT = (acos(-(m2s(0) / lgM2M3)) / M_PI) * 180;
	double aglI = (asin(sin(acos(-(m2s(0,0) / lgM2M3))) * (couplingV / materialV)) / M_PI * 180);

	double xI = -m2s(0);

	double lgei = abs(tan(asin(sin(acos(-(m2s(0) / lgM2M3))) * (couplingV / materialV))) * xI);

	double yI = (lgei / lget) * -m2s(1);
	double zI = (lgei / lget) * -m2s(2);

	MatrixXf m2ts(4,1);
	m2ts << xI, yI, zI, 1;

	MatrixXf m2t;
	m2t = mat3.inverse() * m2ts;

	double r4 = xm1 + couplingH;

	double lgM3M2t = sqrt(pow(xm3 - m2t(0), 2.0) + pow(ym3 - m2t(1), 2.0) + pow(zm3 - m2t(2), 2.0));

	double a3 = (xm3 - m2t(0)) / lgM3M2t;
	double b3 = (ym3 - m2t(1)) / lgM3M2t;
	double c3 = (zm3 - m2t(2)) / lgM3M2t;

	double add11 = pow(xm3, 2.0) + pow(ym3, 2.0);
	double sub11 = add11 - pow(r4, 2.0);
	double mul11 = sub11 * (pow(a3, 2.0) + pow(b3, 2.0)) * 4;

	double discriminant2 = sqrt(pow((xm3 * 2 * a3) + (ym3 * 2 * b3), 2.0) - mul11);
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

	MatrixXf m1PhC(4, 1);
	m1PhC = mat5 * m1;
	MatrixXf m2PhC(4, 1);
	m2PhC = mat5 * m2;
	MatrixXf m3PhC(4, 1);
	m3PhC = mat5 * m3;
	MatrixXf m4PhC(4, 1);
	m4PhC = mat5 * m4;

	double xr3D = abs(m1PhC(0));
	double yr3D = abs(m1PhC(1));
	double zr3D = abs(m1PhC(2));

	double xAlpha3D = abs(m2PhC(0));
	double yAlpha3D = abs(m2PhC(1));
	double zAlpha3D = abs(m2PhC(2));

	double xi3D = abs(m3PhC(0));
	double yi3D = abs(m3PhC(1));
	double zi3D = abs(m3PhC(2));

	double alphaS = acos(xi3D / sqrt(pow(xi3D, 2.0) + pow(yi3D, 2.0) + pow(zi3D, 2.0))) / M_PI * 180;

	vector<double> returnVal {xi3D, yi3D, zi3D, alphaS};

	return returnVal;
}

int main(){
using namespace std;
    auto t = newElipse(337.5, 18);

	cout << t.at(0) << endl;
	cout << t.at(1) << endl;
	cout << t.at(2) << endl;
	cout << t.at(3) << endl;

    return 0;
}