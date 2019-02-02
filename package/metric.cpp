#include "def.h"


//z1就是r，z2就是θ

void metric_KRZ(double spin, double krz_d[], double z1, double z2, double mn[][4]) {
	double r = z1;
	double sqr = r*r;
	double cuber = r*sqr;
	double fourthr = sqr*sqr;

	double th = z2;

	double sinth = sin(th);
	double costh = cos(th);
	double sqsinth = sinth*sinth;
	double sqcosth = costh*costh;

	double sqspin = spin*spin;
	double fourthspin = sqspin*sqspin;
	double r0 = 1 + sqrt(1 - sqspin);
	double sqr0 = r0*r0;
	double cuber0 = sqr0*r0;
	double fourthr0 = sqr0*sqr0;

	double a20 = (2 * sqspin) / cuber0;
	double a21 = -fourthspin / fourthr0 + 0;
	double e0 = (2 - r0) / r0;
	double k00 = sqspin / sqr0;
	double k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
	double w00 = 2 * spin / sqr0;

	double k22 = -sqspin / sqr0;//2017-10-26
	double k23 = sqspin / sqr0;//2017-10-26

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + krz_d[1]*cuber0 / cuber)+(a20*cuber0/cuber+a21*fourthr0/fourthr+k21*cuber0/cuber/(1+k22*(1-r0/r)/(1+k23*(1-r0/r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + krz_d[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26


	mn[0][0] = -(N2 - W*W*sqsinth) / K2;
	mn[0][3] = -W*r*sqsinth;
	mn[1][1] = Sigma*B*B / N2;
	mn[2][2] = Sigma*sqr;
	mn[3][0] = mn[0][3];
	mn[3][3] = K2*sqr*sqsinth;

	return;
}

void metric_KRZ_rderivatives0(double spin, double krz_d[], double z1, double z2, double rdmn[][4])
{

	double r, theta;

	double dr;
	double mnm[4][4];
	double mnp[4][4];

	r = z1;
	theta = z2;
	dr = 0.001*r;


	metric_KRZ(spin, krz_d, r - dr, theta, mnm);
	metric_KRZ(spin, krz_d, r + dr, theta, mnp);

	rdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dr;
	rdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dr;
	rdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dr;
	rdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dr;
	rdmn[3][0] = rdmn[0][3];
	rdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dr;

	return;
}

void metric_KRZ_thderivatives0(double spin, double krz_d[], double z1, double z2, double thdmn[][4])
{
	double r, theta;

	double dtheta;
	double mnm[4][4];
	double mnp[4][4];

	r = z1;
	theta = z2;
	dtheta = 0.01;


	metric_KRZ(spin, krz_d, r, theta - dtheta, mnm);
	metric_KRZ(spin, krz_d, r, theta + dtheta, mnp);

	thdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dtheta;
	thdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dtheta;
	thdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dtheta;
	thdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dtheta;
	thdmn[3][0] = thdmn[0][3];
	thdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dtheta;

	return;
}




void metric_KRZ_rderivatives(double spin, double krz_d[], double z1, double z2, double rdmn[][4]) {

	double r = z1;
	double th = z2;

	double sqr = r*r;
	double cuber = r*sqr;
	double fourthr = sqr*sqr;

	double sinth = sin(th);
	double costh = cos(th);
	double sqsinth = sinth*sinth;
	double sqcosth = costh*costh;

	double sqspin = spin*spin;
	double cubespin = sqspin*spin;
	double fourthspin = sqspin*sqspin;

	double r0 = 1 + sqrt(1 - sqspin);
	double sqr0 = r0*r0;
	double cuber0 = sqr0*r0;
	double fourthr0 = sqr0*sqr0;

	double a20 = (2 * sqspin) / cuber0;
	double a21 = -fourthspin / fourthr0 + 0;
	double e0 = (2 - r0) / r0;
	double k00 = sqspin / sqr0;
	double k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
	double w00 = 2 * spin / sqr0;

	double k22 = -sqspin / sqr0;//2017-10-26
	double k23 = sqspin / sqr0;//2017-10-26

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + krz_d[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + krz_d[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
	//2018-4-4
	//the oringinal N2 was wrong, see the correct KRZ metric

	double rderN2 = (1 - r0 / r)*((2 - r0) / sqr - (2 * (-e0 + k00)*sqr0) / cuber - (3 * krz_d[1]*cuber0) / fourthr) +
		(r0*(1 - (2 - r0) / r + ((-e0 + k00)*sqr0) / sqr + (krz_d[1]*cuber0) / cuber)) / sqr
		+sqcosth*(-3*a20*cuber0/fourthr-4*a21*fourthr0/fourthr/r-3*k21*cuber0/fourthr/ (   1    +      k22*(1-r0/r)/(1+k23*(1-r0/r))  )- 
			k21*cuber0/cuber/  (1+  k22*(1-r0/r)/(1+k23*(1-r0/r ))  ) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))  )  *  ( k22*r0/sqr/(1+  k23*(1-r0/r)  ) - k22*(1-r0/r) / (1 + k23* (1 - r0 / r)) / (1 + k23* (1 - r0 / r)) *k23*r0/sqr  )
			)
		//2018-4-4 modified with the correct metric expression

		/*+
		((-3 * (-0 + fourthspin / fourthr0)*cuber0) / fourthr - (4 * a21*fourthr0) / pow(r, 5))*sqcosth*/;
	double rderB = (-2 * 0 * sqr0) / cuber - (2 * 0 * sqr0*sqcosth) / cuber;
	double rderSigma = (-2 * sqcosth*sqspin) / cuber;
	double rderW = (2 * sqcosth*((krz_d[2] * pow(r0, 3)) / pow(r, 3) + (2 * spin) / pow(r, 2) + (0 * pow(r0, 3)*sqcosth) / pow(r, 3))*sqspin) /
		(pow(r, 3)*pow(1 + (sqcosth*sqspin) / pow(r, 2), 2)) +
		((-3 * krz_d[2] * pow(r0, 3)) / pow(r, 4) - (4 * spin) / pow(r, 3) - (3 * 0 * pow(r0, 3)*sqcosth) / pow(r, 4)) /
		(1 + (sqcosth*sqspin) / pow(r, 2));
	/*	double rderK2 = (2 * cubespin*sqcosth*((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
	(fourthr*pow(1 + (sqcosth*sqspin) / sqr, 2)) +
	(2 * sqcosth*sqspin*((k21*cuber0*sqcosth) / cuber + sqspin / sqr)) /
	(cuber*pow(1 + (sqcosth*sqspin) / sqr, 2)) +
	(spin*((-3 *  0*cuber0) / fourthr - (4 * spin) / cuber - (3 * 0*cuber0*sqcosth) / fourthr)) /
	(r*(1 + (sqcosth*sqspin) / sqr)) - (spin*
	((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
	(sqr*(1 + (sqcosth*sqspin) / sqr)) +
	((-3 * k21*cuber0*sqcosth) / fourthr - (2 * sqspin) / cuber) / (1 + (sqcosth*sqspin) / sqr);*/
	double rderK2 = spin / r*rderW - spin*W / sqr - rderSigma / Sigma / Sigma * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))
		+ (-2 * k00 * sqr0 / cuber - 3 * k21 * cuber0 * sqcosth / fourthr / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))
		- k21 * sqcosth * cuber0 / cuber * (r0 * k22 / sqr / (1 + k23*(1 - r0 / r)) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26

	rdmn[0][0] = (-rderN2 + 2 * rderW*sqsinth*W) / K2 - (rderK2*(-N2 + sqsinth*pow(W, 2))) / pow(K2, 2);

	rdmn[1][1] = (pow(B, 2)*rderSigma) / N2 + (2 * B*rderB*Sigma) / N2 - (pow(B, 2)*rderN2*Sigma) / pow(N2, 2);

	rdmn[2][2] = pow(r, 2)*rderSigma + 2 * r*Sigma;

	rdmn[3][3] = 2 * K2*r*sqsinth + pow(r, 2)*rderK2*sqsinth;

	rdmn[0][3] = -(r*rderW*sqsinth) - sqsinth*W;

	rdmn[3][0] = rdmn[0][3];

	return;
}

void metric_KRZ_thderivatives(double spin, double krz_d[], double z1, double z2, double thdmn[][4]) {
	double r = z1;
	double th = z2;

	double sqr = r*r;
	double cuber = r*sqr;
	double fourthr = sqr*sqr;

	double sinth = sin(th);
	double costh = cos(th); //2018-8-19 之前算cosine居然是开根号。。。不知道以前的raytracing是不是也出过这种问题
	double sqsinth = sinth*sinth;
	double sqcosth = costh*costh;

	double sqspin = spin*spin;
	double cubespin = sqspin*spin;
	double fourthspin = sqspin*sqspin;

	double r0 = 1 + sqrt(1 - sqspin);
	double sqr0 = r0*r0;
	double cuber0 = sqr0*r0;
	double fourthr0 = sqr0*sqr0;

	double a20 = (2 * sqspin) / cuber0;
	double a21 = -fourthspin / fourthr0 + 0;
	double e0 = (2 - r0) / r0;
	double k00 = sqspin / sqr0;
	double k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
	double w00 = 2 * spin / sqr0;

	double k22 = -sqspin / sqr0;//2017-10-26
	double k23 = sqspin / sqr0;//2017-10-26

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + krz_d[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + krz_d[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
//2018-4-4

	double thderN2 = -2 * costh*sinth* (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22* (1 - r0 / r) / (1 + k23 * (1 - r0 / r))));
		//2018-4-4 modified with the correct metric expression
		/* -2 * costh*(((-0 + fourthspin / fourthr0)*cuber0) / cuber + (a21*fourthr0) / fourthr)*sinth */;
	double thderB = (-2 * costh * 0 * sqr0*sinth) / sqr;
	double thderSigma = (-2 * costh*sinth*sqspin) / sqr;
	double thderW = (2 * /*中间好像少了个sqspin???。。。哦没少，写在下一行了*/costh*sinth*((krz_d[2] * cuber0) / cuber + (sqcosth * 0 * cuber0) / cuber + (2 * spin) / sqr)*
		sqspin) / (sqr*pow(1 + (sqcosth*sqspin) / sqr, 2)) -
		(2 * costh * 0 * cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr));
	/*double thderK2 = (2 * costh*cubespin*sinth*((0*cuber0) / cuber + (sqcosth*0*cuber0) / cuber +
	(2 * spin) / sqr)) / (cuber*pow(1 + (sqcosth*sqspin) / sqr, 2)) +
	(2 * costh*sinth*sqspin*((sqcosth*k21*cuber0) / cuber + sqspin / sqr)) /
	(sqr*pow(1 + (sqcosth*sqspin) / sqr, 2)) -
	(2 * costh*k21*cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr)) -
	(2 * costh*0*cuber0*sinth*spin) / (fourthr*(1 + (sqcosth*sqspin) / sqr));*/

	double thderK2 = spin / r * thderW - thderSigma / Sigma / Sigma  * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))
		+ (-2* /*2018-4-8 a factor of 2 was missed...*/costh * sinth *  k21 * cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26


	thdmn[0][0] = (-thderN2 + 2 * pow(sinth, 2)*thderW*W + 2 * costh*sinth*pow(W, 2)) / K2 -
		(thderK2*(-N2 + pow(sinth, 2)*pow(W, 2))) / pow(K2, 2);

	thdmn[1][1] = (2 * B*Sigma*thderB) / N2 - (pow(B, 2)*Sigma*thderN2) / pow(N2, 2) + (pow(B, 2)*thderSigma) / N2;

	thdmn[2][2] = pow(r, 2)*thderSigma;

	thdmn[3][3] = 2 * costh*K2*pow(r, 2)*sinth + pow(r, 2)*pow(sinth, 2)*thderK2;

	thdmn[0][3] = -(r*pow(sinth, 2)*thderW) - 2 * costh*r*sinth*W;

	thdmn[3][0] = thdmn[0][3];
	/*double mnm[4][4];
	double mnp[4][4];
	double dtheta,theta;

	r = z1;
	theta = z2;
	dtheta = 0.01;


	metric_KRZ(spin, krz_d[], r, theta - dtheta, mnm);
	metric_KRZ(spin, krz_d[], r, theta + dtheta, mnp);*/

	return;
}

void metric_KRZ_inverse(double spin, double krz_d[], double z1, double z2, double invg[][4]) {
	double r = z1;
	double th = z2;

	double sqr = r*r;
	double cuber = r*sqr;
	double fourthr = sqr*sqr;

	double sinth = sin(th);
	double costh = sqrt(1 - sinth*sinth);
	double sqsinth = sinth*sinth;
	double sqcosth = costh*costh;

	double sqspin = spin*spin;
	double cubespin = sqspin*spin;
	double fourthspin = sqspin*sqspin;

	double r0 = 1 + sqrt(1 - sqspin);
	double sqr0 = r0*r0;
	double cuber0 = sqr0*r0;
	double fourthr0 = sqr0*sqr0;

	double a20 = (2 * sqspin) / cuber0;
	double a21 = -fourthspin / fourthr0 + 0;
	double e0 = (2 - r0) / r0;
	double k00 = sqspin / sqr0;
	double k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
	double w00 = 2 * spin / sqr0;

	double k22 = -sqspin / sqr0;//2017-10-26
	double k23 = sqspin / sqr0;//2017-10-26
	/*
	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + krz_d[]*cuber0 / cuber) + ((k21 + a20)*cuber0 / cuber + a21*fourthr0 / fourthr)*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26
	*/
	
	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + krz_d[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + krz_d[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26



	invg[0][0] = -(K2 / N2);
	invg[1][1] = N2 / (pow(B, 2)*Sigma);
	invg[2][2] = 1 / (sqr*Sigma);
	invg[3][3] = (-pow(W, 2) + N2 / pow(sin(th), 2)) / (K2*N2*sqr);
	invg[0][3] = -(W / (N2*r));
	invg[3][0] = invg[0][3];


	/*double g[4][4], gg;
	metric_KRZ(spin, krz_d[], z1, z2, g);
	//metric_KRZ_inverse(spin, krz_d[],  w1, w2, invg);
	gg = g[0][0] * g[3][3] - g[0][3] * g[0][3];
	if (std::fabs(gg) > 1e10) {
		invg[0][0] = g[3][3] / gg;
		invg[0][3] = -g[0][3] / gg;
		invg[1][1] = 1 / g[1][1];
		invg[2][2] = 1 / g[2][2];
		invg[3][0] = invg[0][3];
		invg[3][3] = g[0][0] / gg;
	}*/
	return;
}


void find_isco_KRZ(double spin, double krz_d[], double z1, double& isco)
{

	int i, j, casenum = 1, stop = 5, count = 0;
	double detol = 1.0e-5;
	double rll, rul, rinit = z1, rnew, rold;
	double mn[4][4], dmn[4][4];
	double mnp[4][4], mnm[4][4], dmnp[4][4], dmnm[4][4];
	double ep, em, eold, enew, om, omp, omm, omold, omnew;
	double deold, denew;
	double dr = 1.0e-5;
	double sqnu_r;
	double sqspin = spin*spin;

	if (spin>0)
		rll = 1. + sqrt(1. - sqspin);
	else if (spin<0)
		rll = 1. + sqrt(1. - sqspin);
	else
		rll = 1;


	rul = rinit; rold = rul;
	metric_KRZ(spin, krz_d, rold + dr, Pi / 2., mnp);
	metric_KRZ(spin, krz_d, rold - dr, Pi / 2., mnm);
	metric_KRZ_rderivatives(spin, krz_d, rold + dr, Pi / 2., dmnp);
	metric_KRZ_rderivatives(spin, krz_d, rold - dr, Pi / 2., dmnm);
	omp = (-dmnp[0][3] + sqrt(dmnp[0][3] * dmnp[0][3] - dmnp[0][0] * dmnp[3][3])) / dmnp[3][3];
	omm = (-dmnm[0][3] + sqrt(dmnm[0][3] * dmnm[0][3] - dmnm[0][0] * dmnm[3][3])) / dmnm[3][3];
	ep = -(mnp[0][0] + mnp[0][3] * omp) / sqrt(-mnp[0][0] - 2.*mnp[0][3] * omp - mnp[3][3] * omp*omp);
	em = -(mnm[0][0] + mnm[0][3] * omm) / sqrt(-mnm[0][0] - 2.*mnm[0][3] * omm - mnm[3][3] * omm*omm);
	deold = 0.5*(ep - em) / dr;

	do {
		count++;
		if (count>100) {
			printf("No convergence after %i iterations. deold = %.5Le, denew = %.5Le\n", count, deold, denew);
			break;
		}

		rnew = (rll + rul) / 2.;
		metric_KRZ(spin, krz_d, rnew + dr, Pi / 2., mnp);
		metric_KRZ(spin, krz_d, rnew - dr, Pi / 2., mnm);
		metric_KRZ_rderivatives(spin, krz_d, rnew + dr, Pi / 2., dmnp);
		metric_KRZ_rderivatives(spin, krz_d, rnew - dr, Pi / 2., dmnm);
		omp = (-dmnp[0][3] + sqrt(dmnp[0][3] * dmnp[0][3] - dmnp[0][0] * dmnp[3][3])) / dmnp[3][3];
		omm = (-dmnm[0][3] + sqrt(dmnm[0][3] * dmnm[0][3] - dmnm[0][0] * dmnm[3][3])) / dmnm[3][3];
		ep = -(mnp[0][0] + mnp[0][3] * omp) / sqrt(-mnp[0][0] - 2.*mnp[0][3] * omp - mnp[3][3] * omp*omp);
		em = -(mnm[0][0] + mnm[0][3] * omm) / sqrt(-mnm[0][0] - 2.*mnm[0][3] * omm - mnm[3][3] * omm*omm);
		denew = 0.5*(ep - em) / dr;
		//printf("%Le\t%Le\n", denew, rnew);
		if (std::fabs(denew)<std::fabs(detol)) {
			//printf("spin = %Le\tdenew = %Le, deold = %Le\trnew = %Le\n",spin,denew,deold,rnew);
			
			//printf("Normal\n");
			stop = 1;
		}
		else if ((denew*deold)>0.0) {
			if (rnew<rold)
				rul = rnew;
			else if (rnew>rold)
				rll = rnew;
			else
				printf("rold=rnew? rold = %Le, rnew = %Le\n", rold, rnew);
		}
		else if ((denew*deold)<0.0) {
			if (rnew<rold)
				rll = rnew;
			else if (rnew>rold)
				rul = rnew;
			else
				printf("rold=rnew? rold = %Le, rnew = %Le\n", rold, rnew);
		}
		else {
			//printf("Compare enew and eold. eold = %Le, enew = %Le\n", deold, denew);
			break;
		}
		//printf("ep = %Le, em = %Le, rold = %Le, rnew = %Le\n", ep, em, rold, rnew);
		rold = rnew;
		//omold = omnew;
		deold = denew;
	} while (stop == 5);




	/*rold = 6.0;
	rnew = rold;

	do {
	radial_stability_KRZ(spin, krz_d[],   rnew, sqnu_r, om);
	if (sqnu_r>0 && om>0) {
	metric_KRZ(spin, krz_d[],   rnew + dr, Pi / 2., mnp);
	metric_KRZ(spin, krz_d[],   rnew - dr, Pi / 2., mnm);
	metric_KRZ_rderivatives(spin, krz_d[],   rnew + dr, Pi / 2., dmnp);
	metric_KRZ_rderivatives(spin, krz_d[],   rnew - dr, Pi / 2., dmnm);
	omp = (-dmnp[0][3] + sqrt(dmnp[0][3] * dmnp[0][3] - dmnp[0][0] * dmnp[3][3])) / dmnp[3][3];
	omm = (-dmnm[0][3] + sqrt(dmnm[0][3] * dmnm[0][3] - dmnm[0][0] * dmnm[3][3])) / dmnm[3][3];
	ep = -(mnp[0][0] + mnp[0][3] * omp) / sqrt(-mnp[0][0] - 2.*mnp[0][3] * omp - mnp[3][3] * omp*omp);
	em = -(mnm[0][0] + mnm[0][3] * omm) / sqrt(-mnm[0][0] - 2.*mnm[0][3] * omm - mnm[3][3] * omm*omm);
	denew = 0.5*(ep - em) / dr;

	if (std::fabs(denew)<detol) {
	stop = 1;
	break;
	}
	else {
	rold = rnew;
	rnew *= 0.95;
	}
	}
	else {
	stop = 2;
	break;
	}
	} while (stop>2);*/



	if (stop == 1)
		isco = rnew;
	else if (stop == 2) {
		printf("ISCO occurs before minima in E. rold = %Le\n", rold);
		isco = rold;
	}
	else
		isco = 1+sqrt(1-sqspin);

	//printf("stop = %i\n",stop);

	return;
}

void find_isco_KRZ_freq_method(double spin, double krz_d[], double z1, double& isco)
{

	int i, j, casenum = 1, stop = 5, count = 0;
	double detol = 1.0e-5;
	double rll, rul, rinit = z1, rnew, rold;
	double mn[4][4], dmn[4][4];
	double mnp[4][4], mnm[4][4], dmnp[4][4], dmnm[4][4];
	double ep, em, eold, enew, om, omp, omm, omold, omnew;
	double deold, denew;
	double dr = 1.0e-5;
	double sqnu_r;
	double sqnu_th;
	double sqspin = spin*spin;

	if (spin>0)
		rll = 1. + sqrt(1. - sqspin);
	else if (spin<0)
		rll = 1. + sqrt(1. - sqspin);
	else
		rll = 1;


	rold = 12.0;//It used to start from 6.0
	rnew = rold;

	//------------		SANDBOX		--------------//
	for (rnew = rold; rnew>(rll*1.001); rnew -= 0.001){
		metric_KRZ(spin, krz_d, rnew, Pi / 2., mnp);
		metric_KRZ_rderivatives(spin, krz_d, rnew, Pi / 2., dmnp);
		omp = (-dmnp[0][3] + sqrt(dmnp[0][3] * dmnp[0][3] - dmnp[0][0] * dmnp[3][3])) / dmnp[3][3];
		ep = -(mnp[0][0] + mnp[0][3] * omp) / sqrt(-mnp[0][0] - 2.*mnp[0][3] * omp - mnp[3][3] * omp*omp);
		printf("%Le %Le\n", rnew, ep);
	}
	//abort();


	do {
		radial_stability_KRZ(spin, krz_d, rnew, sqnu_r, om);
		vertical_stability_KRZ(spin, krz_d, rnew, sqnu_th, om);
		//printf("sqnu_r = %Le\tr = %Le\n", sqnu_r, rnew);
		if (sqnu_r > 0 && sqnu_th > 0 && om > 0) {
			metric_KRZ(spin, krz_d, rnew + dr, Pi / 2., mnp);
			metric_KRZ(spin, krz_d, rnew - dr, Pi / 2., mnm);
			metric_KRZ_rderivatives(spin, krz_d, rnew + dr, Pi / 2., dmnp);
			metric_KRZ_rderivatives(spin, krz_d, rnew - dr, Pi / 2., dmnm);
			omp = (-dmnp[0][3] + sqrt(dmnp[0][3] * dmnp[0][3] - dmnp[0][0] * dmnp[3][3])) / dmnp[3][3];
			omm = (-dmnm[0][3] + sqrt(dmnm[0][3] * dmnm[0][3] - dmnm[0][0] * dmnm[3][3])) / dmnm[3][3];
			ep = -(mnp[0][0] + mnp[0][3] * omp) / sqrt(-mnp[0][0] - 2.*mnp[0][3] * omp - mnp[3][3] * omp*omp);
			em = -(mnm[0][0] + mnm[0][3] * omm) / sqrt(-mnm[0][0] - 2.*mnm[0][3] * omm - mnm[3][3] * omm*omm);
			denew = 0.5*(ep - em) / dr;

			if (std::fabs(denew) < detol) {
				stop = 1;
				break;
			}
			else {
				//printf("ep = %Le, em = %Le, rold = %Le, rnew = %Le\n", ep, em, rold, rnew);
				rold = rnew;
				rnew *= 0.999;//It was rnew*=0.95 before.
			}
		}
		else {
			stop = 2;
			break;
		}
	} while (stop>2);



	if (stop == 1)
		isco = rnew;
	else if (stop == 2) {
		printf("ISCO occurs before minima in E. rold = %Le sqnu_r = %Le sqnu_th = %Le omega = %Le\n", rold, sqnu_r, sqnu_th, om);
		isco = rold;
	}
	else
		isco = 0;

	//printf("stop = %i\n",stop);

	return;
}


void radial_stability_KRZ(double spin, double krz_d[], double z1, double& sqfreq_r, double& omega)
{
	int i;
	double r, dr, th, dth;
	double en, elle;
	double mn[4][4], dmn[4][4];
	double tdot, Veff[3];

	char filename_o[128];
	FILE *fout;

	r = z1;
	dr = 0.01*r;
	//dth = 0.01;

	//sprintf(filename_o,"vertistab_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
	//fout = fopen(filename_o,"w");

	//do{
	th = Pi / 2;
	metric_KRZ(spin, krz_d, r, th, mn);
	metric_KRZ_rderivatives(spin, krz_d, r, th, dmn);
	omega = (-dmn[0][3] + sqrt(dmn[0][3] * dmn[0][3] - dmn[0][0] * dmn[3][3])) / dmn[3][3];
	en = -(mn[0][0] + mn[0][3] * omega) / sqrt(-mn[0][0] - 2.*mn[0][3] * omega - mn[3][3] * omega*omega);
	elle = (mn[0][3] + mn[3][3] * omega) / sqrt(-mn[0][0] - 2.*mn[0][3] * omega - mn[3][3] * omega*omega);

	tdot = (en*mn[3][3] + elle*mn[0][3]) / (mn[0][3] * mn[0][3] - mn[0][0] * mn[3][3]);

	for (i = 0; i<3; i++) {
		r = z1 - dr + (double)i*dr;
		metric_KRZ(spin, krz_d, r, th, mn);
		Veff[i] = (en*en*mn[3][3] + 2 * en*elle*mn[0][3] + elle*elle*mn[0][0]) / (mn[0][3] * mn[0][3] - mn[0][0] * mn[3][3]);
	}

	sqfreq_r = -0.5*(Veff[2] + Veff[0] - 2 * Veff[1]) / (dr*dr*mn[1][1] * tdot*tdot);

	//fprintf(fout,"%.10Le %.10le\n",r,sqrt(sqfreq_th));

	//r-=dr;

	//}while(sqfreq_th>0 && r>isco);


	//fclose(fout);
	return;

}

void vertical_stability_KRZ(double spin, double krz_d[], double z1, double& sqfreq_th, double omega) {
	int i;
	double r, dth, th;
	double en, elle;
	double mn[4][4], dmn[4][4];
	double tdot, Veff[3];


	//char filename_o[128];
	//FILE *fout;

	r = z1;
	th = Pi / 2;
	dth = 0.01*th;
	//dth = 0.01;

	//sprintf(filename_o,"vertistab_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
	//fout = fopen(filename_o,"w");

	//do{

	metric_KRZ(spin, krz_d, r, th, mn);
	metric_KRZ_thderivatives(spin, krz_d, r, th, dmn);

	en = -(mn[0][0] + mn[0][3] * omega) / sqrt(-mn[0][0] - 2.*mn[0][3] * omega - mn[3][3] * omega*omega);
	elle = (mn[0][3] + mn[3][3] * omega) / sqrt(-mn[0][0] - 2.*mn[0][3] * omega - mn[3][3] * omega*omega);

	tdot = (en*mn[3][3] + elle*mn[0][3]) / (mn[0][3] * mn[0][3] - mn[0][0] * mn[3][3]);

	for (i = 0; i < 3; i++) {
		th = Pi / 2 - dth + (double)i*dth;
		metric_KRZ(spin, krz_d, r, th, mn);
		Veff[i] = (en*en*mn[3][3] + 2 * en*elle*mn[0][3] + elle*elle*mn[0][0]) / (mn[0][3] * mn[0][3] - mn[0][0] * mn[3][3]);
	}

	sqfreq_th = -0.5*(Veff[2] + Veff[0] - 2 * Veff[1]) / (dth*dth*mn[2][2] * tdot*tdot);

	//fprintf(fout,"%.10Le %.10le\n",r,sqrt(sqfreq_th));

	//r-=dr;

	//}while(sqfreq_th>0 && r>isco);


	//fclose(fout);
	return;
}
