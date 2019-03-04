#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <fstream>
using namespace std;

const double eps = 1e-15;//double是16位小数
const double myinf = 1e20;

const double Pi = 3.141592653589793238462643383279502884197;
const double Piby2 = 1.57079632679489655799898173427209258;


void metric(double spin, double defpar[], double z1, double z2, double mn[][4]);
void metric_rderivatives(double spin, double defpar[], double z1, double z2, double rdmn[][4]);
void metric_thderivatives(double spin, double defpar[], double z1, double z2, double thdmn[][4]);
void metric_inverse(double spin, double defpar[], double z1, double z2, double invg[][4]);
void Christoffel(double spin, double defpar[], double w1, double w2, double CS[][4][4]);


#endif