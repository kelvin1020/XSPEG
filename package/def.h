#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <fstream>
using namespace std;

const double Pi  = 3.141592653589793;
const double Piby2 = 1.5707963267948966192;

//double mymax(double a, double b);
//void diffeqs(double b, double vars[], double diffs[],double spin, double krz_d[]);
//void intersection(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double x_d[]);
//void raycovg(double spin, double epsi3, double a13, double a22, double a52, double robs, double pobs, double iobs, double xin, double germ, double traces[]);

/****     KRZ metric    ****/

void metric_KRZ(double spin, double krz_d[], double z1, double z2, double mn[][4]);
void metric_KRZ_rderivatives(double spin, double krz_d[], double z1, double z2, double rdmn[][4]);
void metric_KRZ_thderivatives(double spin, double krz_d[], double z1, double z2, double thdmn[][4]);
void metric_KRZ_inverse(double spin, double krz_d[], double z1, double z2, double invg[][4]);
void Christoffel_KRZ(double spin, double krz_d[], double w1, double w2, double CS[][4][4]);

void radial_stability_KRZ(double spin, double krz_d[], double z1, double& sqfreq_r, double& omega);
void find_isco_KRZ(double spin, double krz_d[],  double z1, double& isco);
void find_isco_KRZ_freq_method(double spin, double krz_d[], double z1, double& isco);
//void jacobian_KRZ(double spin, double krz_d[],  double robscur, double pobs, double iobs, double isco, double xin, double robs, double& jac);
//void redshift_KRZ(double spin, double krz_d[],  double radius, double ktkp, double& gg);

//void raytrace_KRZ(double spin, double krz_d[], double xobs, double yobs, double iobs, double xin, double traces[5]);
//void rayprecise_KRZ(double spin, double krz_d[],  double robs, double germtol, double iobs, double xin, double isco, double pobs, double traced[5]);

//void xyfromrphi(double robscur, double pobs, double iobs, double & xobs, double & yobs);
/*new functions*/
//double gradient(double spin, double krz_d[], double robs, double pobs, double iobs, double xin, double robscur, double fmiss);
void vertical_stability_KRZ(double spin, double krz_d[],  double z1, double& sqfreq_th, double omega);
#endif