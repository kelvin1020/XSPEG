#ifndef _DEF_H
#include "def.h"
#endif


void metric(double spin /* the black hole spin */, double defpar[]/* deformation parameter in the custom metric */,
	double z1/* the r component of Boyer Lindquist coordinates */, double z2/* the theta component of Boyer Lindquist coordinates */, 
	double mn[][4]/* the metric elements */) {

	/* write your metric here, below is an example of KRZ metric */
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

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + defpar[1]*cuber0 / cuber)+(a20*cuber0/cuber+a21*fourthr0/fourthr+k21*cuber0/cuber/(1+k22*(1-r0/r)/(1+k23*(1-r0/r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + defpar[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26


	mn[0][0] = -(N2 - W*W*sqsinth) / K2;
	mn[0][3] = -W*r*sqsinth;
	mn[1][1] = Sigma*B*B / N2;
	mn[2][2] = Sigma*sqr;
	mn[3][0] = mn[0][3];
	mn[3][3] = K2*sqr*sqsinth;

	return;
}



void metric_rderivatives(double spin /* the black hole spin */, double defpar[]/* deformation parameter in the custom metric */,
	double z1/* the r component of Boyer Lindquist coordinates */, double z2/* the theta component of Boyer Lindquist coordinates */,
	double rdmn[][4]/* metric derivatives with respect to r */) {
	
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

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + defpar[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + defpar[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;

	double rderN2 = (1 - r0 / r)*((2 - r0) / sqr - (2 * (-e0 + k00)*sqr0) / cuber - (3 * defpar[1] * cuber0) / fourthr) +
		(r0*(1 - (2 - r0) / r + ((-e0 + k00)*sqr0) / sqr + (defpar[1] * cuber0) / cuber)) / sqr
		+ sqcosth*(-3 * a20*cuber0 / fourthr - 4 * a21*fourthr0 / fourthr / r - 3 * k21*cuber0 / fourthr / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))) -
			k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))  *  (k22*r0 / sqr / (1 + k23*(1 - r0 / r)) - k22*(1 - r0 / r) / (1 + k23* (1 - r0 / r)) / (1 + k23* (1 - r0 / r)) *k23*r0 / sqr)
			);

	double rderB = (-2 * 0 * sqr0) / cuber - (2 * 0 * sqr0*sqcosth) / cuber;
	double rderSigma = (-2 * sqcosth*sqspin) / cuber;
	double rderW = (2 * sqcosth*((defpar[2] * pow(r0, 3)) / pow(r, 3) + (2 * spin) / pow(r, 2) + (0 * pow(r0, 3)*sqcosth) / pow(r, 3))*sqspin) /
		(pow(r, 3)*pow(1 + (sqcosth*sqspin) / pow(r, 2), 2)) +
		((-3 * defpar[2] * pow(r0, 3)) / pow(r, 4) - (4 * spin) / pow(r, 3) - (3 * 0 * pow(r0, 3)*sqcosth) / pow(r, 4)) /
		(1 + (sqcosth*sqspin) / pow(r, 2));

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

void metric_thderivatives(double spin /* the black hole spin */, double defpar[]/* deformation parameter in the custom metric */,
	double z1/* the r component of Boyer Lindquist coordinates */, double z2/* the theta component of Boyer Lindquist coordinates */,
	double thdmn[][4]/* metric derivatives with respect to theta */) {
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

	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + defpar[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + defpar[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
//2018-4-4

	double thderN2 = -2 * costh*sinth* (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22* (1 - r0 / r) / (1 + k23 * (1 - r0 / r))));
	double thderB = (-2 * costh * 0 * sqr0*sinth) / sqr;
	double thderSigma = (-2 * costh*sinth*sqspin) / sqr;
	double thderW = (2 *costh*sinth*((defpar[2] * cuber0) / cuber + (sqcosth * 0 * cuber0) / cuber + (2 * spin) / sqr)*
		sqspin) / (sqr*pow(1 + (sqcosth*sqspin) / sqr, 2)) -
		(2 * costh * 0 * cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr));
	

	double thderK2 = spin / r * thderW - thderSigma / Sigma / Sigma  * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))
		+ (-2* costh * sinth *  k21 * cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26


	thdmn[0][0] = (-thderN2 + 2 * pow(sinth, 2)*thderW*W + 2 * costh*sinth*pow(W, 2)) / K2 -
		(thderK2*(-N2 + pow(sinth, 2)*pow(W, 2))) / pow(K2, 2);

	thdmn[1][1] = (2 * B*Sigma*thderB) / N2 - (pow(B, 2)*Sigma*thderN2) / pow(N2, 2) + (pow(B, 2)*thderSigma) / N2;

	thdmn[2][2] = pow(r, 2)*thderSigma;

	thdmn[3][3] = 2 * costh*K2*pow(r, 2)*sinth + pow(r, 2)*pow(sinth, 2)*thderK2;

	thdmn[0][3] = -(r*pow(sinth, 2)*thderW) - 2 * costh*r*sinth*W;

	thdmn[3][0] = thdmn[0][3];

	return;
}

void metric_inverse(double spin, double defpar[], double z1, double z2, double invg[][4]) {
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
	double N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + defpar[1]*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
	double B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
	double Sigma = 1 + sqspin*sqcosth / sqr;
	double W = (w00*sqr0 / sqr + defpar[2] * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
	double K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;//2017-10-26



	invg[0][0] = -(K2 / N2);
	invg[1][1] = N2 / (pow(B, 2)*Sigma);
	invg[2][2] = 1 / (sqr*Sigma);
	invg[3][3] = (-pow(W, 2) + N2 / pow(sin(th), 2)) / (K2*N2*sqr);
	invg[0][3] = -(W / (N2*r));
	invg[3][0] = invg[0][3];

	return;
}
