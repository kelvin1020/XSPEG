#include "def.h"
#include <stdlib.h>
#include <iostream>
#include <string>
//#include <time.h>
using namespace std;
#ifndef max
#define max(a,b) (( (a) >= (b)) ? (a) : (b))
#endif
#define N 8 //8¸ö±äÁ¿

void equations(double var[],double diff[],double spin,double krz_d[]);



int main(int argc, char *argv[])
{//RK45µÄ²ÎÊý
	double a1 = 1.0 / 4.0;
	double b1 = 3.0 / 32.0;
	double b2 = 9.0 / 32.0;
	double c1 = 1932.0 / 2197.0;
	double c2 = -7200.0 / 2197.0;
	double c3 = 7296.0 / 2197.0;
	double d1 = 439.0 / 216.0;
	double d2 = -8.0;
	double d3 = 3680.0 / 513.0;
	double d4 = -845.0 / 4104.0;
	double e1 = -8.0 / 27.0;
	double e2 = 2.0;
	double e3 = -3544.0 / 2565.0;
	double e4 = 1859.0 / 4104.0;
	double e5 = -11.0 / 40.0;
	double x1 = 25.0 / 216.0;
	double x2 = 0.0;
	double x3 = 1408.0 / 2565.0;
	double x4 = 2197.0 / 4104.0;
	double x5 = -1.0 / 5.0;
	double z1 = 16.0 / 135.0;
	double z2 = 0.0;
	double z3 = 6656.0 / 12825.0;
	double z4 = 28561.0 / 56430.0;
	double z5 = -9.0 / 50.0;
	double z6 = 2.0 / 55.0;

	double spin, spin2, krz_d[10] = { 0 };
	double isco, xin;

	double robs_i, robs_f;
	double pstep;



	int i, j, k, m;
	int ii,jj;
	char filename_o[128];

	FILE *foutput, *finput;

	double r, th, t, phi, ur, ut, uth, uphi, tau, dtau; //×ø±ê£¬4ËÙ¶È£¬¹ÌÓÐÊ±

	//double rnew, thnew, tnew, phinew, urnew, utnew, uthnew, uphinew, dtau;

	double E, Lz;

	double f_t, f_phi, f_theta, F_theta, F_r,F_t,F_phi;
	double F_thetanew, F_rnew;

	double g[4][4] = { 0 }, gnew[4][4];
	double Gamma[4][4][4], Gammanew[4][4][4];
	double u[4], unew[4];
	double k1[8], k2[8], k3[8], k4[8],k5[8],k6[8];//8¸ö±äÁ¿µÄË³ÐòÊÇt,r,\theta,\phi,u^r,u^\theta, u^t, u^\phi
	double var[8];
	double ita,mu;//¼ÆËã¹ý³ÌÖÐËÄËÙ¶ÈµÄÄ£\ita£» ÓÐÖÊÁ¿£ºmu=-1£¬ÎÞÖÊÁ¿£º mu=0

	/***************************¸Ä³ÉÃë¼ä¸ôÐèÒª¸ÄµÄ²¿·Ö********************/

	double dt;
	double Grav, clight, Msol, M, dt_sec = 0.1;
	double t_sec = 0;
	
	Grav = 6.674e-11;//#ÒýÁ¦³£Êý
	clight = 2.998e8;//#¹âËÙ
	Msol = 1.989e30;  //#Ì«ÑôÖÊÁ¿£¬ÒÔÇ§¿Ë×öµ¥Î»
//input argument: 1. M, 2. spin, 3. E, 4. Lz, 5. Q, 6. r0, 7. tottime, 8-10. krz_d
	M = atof(argv[1]);
	//dt_sec = dt*M*Msol*Grav / clight / clight / clight

	dt = dt_sec*clight*clight*clight / M / Msol / Grav;

	//dt = 1;

	/***************************¸Ä³ÉÃë¼ä¸ôÐèÒª¸ÄµÄ²¿·Ö********************/

	mu = -1;
	tau = 0;
	dtau = dt;//²½³¤
	double tottime=atof(argv[7]);
	spin = atof(argv[2]);//+0.001*0.5;
	//for (spin = 0.55428452790227312;spin <  0.55428452790227312 +0.01;spin += 0.05) {
		printf("Metric parameters: M=%f, spin=%f \nKRZ metric\n",M, spin);
		krz_d[1] = atof(argv[8]),krz_d[2]=atof(argv[9]), krz_d[3]=atof(argv[10]);//+0.001*0.2;
//		for (krz_d[2] = 0.0; krz_d[2] <= 0.0; krz_d[2] += 0.1) {
			printf("delta: %f, %f, %f\n",krz_d[1], krz_d[2],krz_d[3]);
			t = 0;
			//r = 13.0;
			th = Piby2;
			phi = 0;


			//ut = ( E*g[3][3]+Lz*g[0][3] ) / ( g[0][3]*g[0][3] - g[0][0]*g[3][3] );
			//uphi = (E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]) ;
			//uth = sqrt((-1 - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));

			/********************¡ýÔ²¹ìµÀµÄÄÜÁ¿ºÍ½Ç¶¯Á¿ from https://arxiv.org/pdf/1105.2959.pdf £¨ºÃÏñÊÇ´íµÄ¡£¡££©¡ý******************/
			/*	Lz = fabs(spin*spin + 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) - 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, corotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin - sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) - 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, corotating

			Lz = fabs(spin*spin - 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) + 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, counterrotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin + sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) + 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, counterrotating
			*/
			/********************¡üfrom https://arxiv.org/pdf/1105.2959.pdf £¨ºÃÏñÊÇ´íµÄ¡£¡£¡££©¡ü******************/

			/********************¡ýÔ²¹ìµÀµÄÄÜÁ¿ºÍ½Ç¶¯Á¿ from http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf £¨Õâ¸öÊÇ¶ÔµÄ£©¡ý******************/

			//Lz = ( sqrt(r) -2*spin/r + spin*spin/sqrt(r*r*r) ) / sqrt(1 - 3 / r + 2 * spin / sqrt(r*r*r));
			//E = (1-2/r + spin /sqrt(r*r*r) ) / sqrt(1-3/r + 2*spin/sqrt(r*r*r) );

			/********************¡üfrom http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf ¡ü******************/


			double horizon, r0, th0 = Piby2, phi0 = 0, t0 = 0;//³õÊ¼Ìõ¼þ
			double ecc, p;//eccentricityºÍrectum 
			double rmax, rmin, invgmax[4][4], invgmin[4][4];//rµÄÉÏÏÂÏÞÒÔ¼°¸Ã´¦µÄg
			double invg[4][4];
			double EoverL;//ÓÉeºÍpËãEºÍLzµÄÖÐ¼ä±äÁ¿
			double EoverL2, E2, L2;
			//double iota = Pi / 6;
			double Q,Q0;//carter constant
			//horizon = 1 + sqrt(1 - spin*spin);
			/***************************³àµÀÃæÉÏÓÉe,p¾ö¶¨µÄ¹ìµÀ¡ý******************/
			//find_isco_KRZ(spin, krz_d, 100., isco);
			horizon= 1 +  sqrt(1-spin*spin);//atof(argv[11]);
			//for (ecc = 0.41149551640029858;ecc <= 0.41149551640029858;ecc = ecc + 0.004) {
			/*for (ecc = 0.5;ecc <= 0.5;ecc = ecc + 0.004) {
				//	for (p = 6.4825501282607396;p <= 6.4825501282607396;p = p + 0.04) {
				for (p = 6;p <= 6;p = p + 0.04) {*/
					
			//ecc = 0.43173473149300562,p= 7.1717260434110983;
					E = atof(argv[3]), Lz=atof(argv[4]), Q=atof(argv[5]);//python½â³öÀ´µÄ   and note that it's not useful actually, Q is calculated below, and for non-Kerr case this is just initial Q
					th0 = Piby2;//³àµÀÃæÉÏcarter constantºÍtheta·½ÏòËÙ¶È¹ØÏµ±È½Ï¼òµ¥
					r0 = atof(argv[6]);//¿´¿´È¡³õÊ¼ÔÚrmax»áÔõÃ´Ñù£¿
					Q0 = Q;
			/*
					rmax = p / (1 - ecc);
					rmin = p / (1 + ecc);
					r0 = rmax;
					metric_KRZ_inverse(spin, krz_d, rmax, th0, invgmax);
					metric_KRZ_inverse(spin, krz_d, rmin, th0, invgmin);


					EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / (invgmax[0][0] - invgmin[0][0]);
					Lz = sqrt((invgmax[3][0] - invgmin[3][0]) / (EoverL*EoverL*(invgmin[3][0] * invgmax[0][0] - invgmax[3][0] * invgmin[0][0]) + (invgmin[3][0] * invgmax[3][3] - invgmax[3][0] * invgmin[3][3])));

					E = EoverL*Lz;

					*/
					//EoverL2 = ((invgmax[3][0] - invgmin[3][0]) - sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / (invgmax[0][0] - invgmin[0][0]);
					//L2 = sqrt((invgmax[3][0] - invgmin[3][0]) / (EoverL*EoverL*(invgmin[3][0] * invgmax[0][0] - invgmax[3][0] * invgmin[0][0]) + (invgmin[3][0] * invgmax[3][3] - invgmax[3][0] * invgmin[3][3])));
					//
					//E2 = abs(EoverL2*Lz);
					//if (E2 < E) {
					//	E = E2;
					//}
					r = r0;
					th = th0;
					t = t0;
					phi = phi0;
					metric_KRZ(spin, krz_d, r, th, g);
					metric_KRZ_inverse(spin, krz_d, r, th, invg);



					//printf("E=%.6f Lz=%.6f Q=%.6f\n",  E, Lz, Q);
					/********************¡ý¶þ·ÖÕÒE,Lz¡ý******************/
					/*
					if (abs(ecc - 0.0) < 1e-6) {
						r = p;r0 = p;
						double upE = 1, downE = 0.9, eps = 1e-10;
						double invg[4][4];
						double curf, upf, downf;
						Christoffel_KRZ(spin, krz_d, r, th, Gamma);
						metric_KRZ_inverse(spin, krz_d, r, th, invg);

						E = upE;
						Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

						ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						upf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;


						E = downE;
						Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

						ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
						downf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;


						for (;;) {
							E = 0.5*(upE + downE);
							Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];

							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							curf = Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi;

							if (curf*downf < 0) {
								upE = E;
								upf = curf;
							}
							else {
								downE = E;
								downf = curf;
							}
							if (abs(upE - downE) < eps) {
								E = 0.5*(upE + downE);
								Lz = (invg[0][3] * E + sqrt((invg[0][3] * E)*(invg[0][3] * E) - invg[3][3] * (invg[0][0] * E*E + 1))) / invg[3][3];
								break;
							}

						}

					}*/

					/********************¡ü¶þ·ÖÕÒÔ²¹ìµÀ ¡ü******************/
					/*
					double testE, testL, testut, testup, utoverup;
					if (abs(ecc - 0.0) < 1e-6) {
						///½â·½³ÌÕÒÔ²¹ìµÀ¡ý
						//p = 10;
						Christoffel_KRZ(spin, krz_d, p, th, Gamma);
						metric_KRZ(spin, krz_d, p, th, g);
						utoverup = (-Gamma[1][0][3] + sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
						testup = sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
						testut = testup*utoverup;
						testE = -g[0][0] * testut - g[0][3] * testup;
						testL = g[0][3] * testut + g[3][3] * testup;
						E = testE;
						Lz = testL;
						//   ½â·½³ÌÕÒÔ²¹ìµÀ¡ü
					}
					*/

					/*looking for circular orbit¡ý*/
					/*
					metric_KRZ(spin, krz_d, r, th, g);
					Christoffel_KRZ(spin, krz_d, r, th, Gamma);
					sprintf(filename_o, "test.dat");
					foutput = fopen(filename_o, "w");

					for (E = 0.5;E < 1;E += 0.001) {
						for (Lz = 0;Lz < 10;Lz += 0.01) {
							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							fprintf(foutput, "%.6f \t ", Gamma[1][0][0] * ut*ut + 2 * Gamma[1][0][3] * ut*uphi + Gamma[1][3][3] * uphi*uphi);

						}
						for (Lz = 0;Lz < 10;Lz += 0.01) {
							ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
							fprintf(foutput, "%.6f \t ", g[0][0] * ut*ut + 2 *g[0][3] * ut*uphi + g[3][3] * uphi*uphi);

						}
						fprintf(foutput, "\n");
						printf("%.2f %.2f \n", E, Lz);

					}
					fclose(foutput);
					abort();

					/*looking for circular orbit¡ü*/

					ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					//uth = sqrt(Q) / g[2][2];
					//ur = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[2][2] * uth*uth) / (g[1][1]));
					ur = 0;
					uth = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));
					//uth = 0;
					//ur = 0;
					Q = uth*uth*g[2][2] * g[2][2];
					E = -g[0][0] * ut - g[0][3] * uphi;
					Lz = g[0][3] * ut + g[3][3] * uphi;
					u[0] = ut;
					u[1] = ur;
					u[2] = uth;
					u[3] = uphi;

					//printf("template E=%.6f Lz=%.6f\n", ecc, p, E, Lz);

					//sprintf(filename_o, "circular_trace_spin%.2f_d%.2f_r%.2f.dat", spin, krz_d,r);

					printf("\nE=%.6f Lz=%.6f Q=%.6f\n", E, Lz, Q);
					printf("starting condition:\n t=%f, r=%f, th=%f, phi=%f\n ut=%f, ur=%f, uth=%f, uphi=%f\n\n",t,r,th,phi,ut,ur,uth,uphi);
					/***************************¸Ä³ÉÃë¼ä¸ôÐèÒª¸ÄµÄ²¿·Ö********************/

					sprintf(filename_o,"ORBCAR");
					//sprintf(filename_o, "trace_M%.0f_spin%.6f_d1%.6f_d2%.6f_d3%.6f.dat", M, spin, krz_d[1], krz_d[2], krz_d[3]);


					/***************************¸Ä³ÉÃë¼ä¸ôÐèÒª¸ÄµÄ²¿·Ö********************/

					foutput = fopen(filename_o, "w");
					//Zprintf("index\t tau\t t\t r\t theta\t phi\t ut\t ur\t uth\t uphi\t ita\t E\t Lz\n");

					var[0] = t;	var[1] = r;	var[2] = th;	var[3] = phi;	var[4] = ur;	var[5] = uth;	var[6] = ut;	var[7] = uphi;
					double E0 = E, L0 = Lz;
					double h;//RK45, adaptive step
					double diff[8];
					double vars_temp[8];//temp used in RK45
					double z[8], y[8];//used to estimate error
					double err, maxerr=1e-10, minerr=1e-11;//¼ÆËã¹ý³ÌÖÐµÄÎó²î£¬ÈÝÐíµÄ×î´óÎó²î£¬ÈÝÐíµÄ×îÐ¡Îó²î
					int check=0;//±ê¼ÇÓÐÃ»ÓÐ³¬¹ý×î´óÎó²î»òÐ¡ÓÚ×îÐ¡Îó²î
					int index = 0;
					h = dtau;
					tau = 0;
					
					double percentage=0;//monitoring the progress
					//time_t optime = time(NULL);
					for(;t_sec<tottime;){
						if(t_sec/tottime*100>percentage+1.0){
							percentage=t_sec/tottime*100;
                                                        printf("\033[1A"); // move cursor one line up
                                                        printf("\033[K");   // delete till end of line

							printf("Orbit computed %.0f\%  \n", percentage);
						}
						if(var[1]<horizon) break;
						t_sec = var[0] * M*Msol*Grav / clight / clight / clight;

						metric_KRZ(spin, krz_d, var[1], var[2], g);
						uth = var[5];ut = var[6];uphi = var[7];ur = var[4];
//						system("pause");
						ita = g[0][0] * ut*ut + g[1][1] * ur*ur + g[2][2] * uth*uth + g[3][3] * uphi*uphi + 2 * g[3][0] * uphi*ut;
						Q = uth*g[2][2] * uth*g[2][2] + cos(var[2])*cos(var[2])*(spin*spin*(ita*ita - E*E) + Lz*Lz / sin(var[2]) / sin(var[2]));
						E = -g[0][0] * ut - g[0][3] * uphi;
						Lz = g[0][3] * ut + g[3][3] * uphi;
						//printf("%f\n",Lz);getchar();
						//printf("%.10f\n",Lz);
						//printf("%lf\n",Lz);
						//printf("%f\n",Q);getchar();
                                                //printf("%.10f\n",Q);
                                                //printf("%lf\n",Q);


						check = 0;
						equations(var, diff,spin,krz_d);
						for (i = 0; i < N; i++)
						{
							k1[i] = h*diff[i];
							vars_temp[i] = var[i] + a1*k1[i];
						}

						equations(vars_temp, diff,spin,krz_d);
						for (i = 0; i < N; i++)
						{
							k2[i] = h*diff[i];
							vars_temp[i] = var[i] + b1*k1[i] + b2*k2[i];
						}

						equations(vars_temp, diff,spin,krz_d);
						for (i = 0; i < N; i++)
						{
							k3[i] = h*diff[i];
							vars_temp[i] = var[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
						}

						equations(vars_temp, diff,spin,krz_d);
						for (i = 0; i < N; i++)
						{
							k4[i] = h*diff[i];
							vars_temp[i] = var[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
						}

						equations(vars_temp, diff,spin,krz_d);
						for (i = 0; i < N; i++)
						{
							k5[i] = h*diff[i];
							vars_temp[i] = var[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
						}

						equations(vars_temp, diff,spin,krz_d);
						for (i = 0; i < N; i++)
							k6[i] = h*diff[i];

						for (i = 0; i < N; i++)
						{
							y[i] = var[i] + x1*k1[i] + x2*k2[i] + x3*k3[i] + x4*k4[i] + x5*k5[i];
							z[i] = var[i] + z1*k1[i] + z2*k2[i] + z3*k3[i] + z4*k4[i] + z5*k5[i] + z6*k6[i];
							err = fabs((y[i] - z[i]) / max(fabs(var[i]), fabs(y[i]) ) );
							if (err > maxerr) {
								check = 1;//ÓÐÐ©Îó²îÌ«´ó
							}
							else if (err < minerr&&check != 1) {
								check = -1;//ËùÓÐÎó²î¶¼Ì«Ð¡
							}
						}

						if (check == 1) {
							h /= 1.1;
						}
						else if (check == -1) {
							if (check != 1) {
								fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[6], var[4], var[5], var[7], /*F_t*/(z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6]), /*F_r*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]), /*F_theta*/z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5], /*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
								//printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[6]/*ut*/, var[4]/*ur*/, var[5]/*uth*/, var[7]/*uphi*/, ita - mu, (E - E0) / E0, (Lz - L0) / L0, Q);
								fflush(foutput);
							}
							tau = tau + h;
							h *= 1.1;
							index++;
							for (i = 0; i < N; i++)
								var[i] = y[i];
						}
						else {
							if (check != 1) {
								fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[6], var[4], var[5], var[7], /*F_t*/(z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6]), /*F_r*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]), /*F_theta*/z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5], /*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
								//printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n", index, t_sec, tau, var[0], var[1], var[2], var[3], var[6]/*ut*/, var[4]/*ur*/, var[5]/*uth*/, var[7]/*uphi*/, ita - mu, (E - E0) / E0, (Lz - L0) / L0, Q);
								fflush(foutput);
							}
							tau = tau + h;
							index++;
							for (i = 0; i < N; i++)
								var[i] = y[i];
						}
						

					}
					fclose(foutput);
					printf("\033[1A"); // move cursor one line up
                                        printf("\033[K");   // delete till end of line
                                        printf("Orbit computed 100\%  \n");

					//time_t edtime = time(NULL);
					//cout << (edtime - optime) << endl;
					//system("pause");
				
			
		//}

	//}
	return 0;

}

void equations(double var[], double diff[],double spin,double krz_d[]) {//ÓÉ8¸ö±äÁ¿ËãÐ±ÂÊ,×¢Òâ±äÁ¿µÄË³ÐòÊÇ£º0:t, 1:r, 2:th, 3:phi, 4:ur, 5:uth, 6:ut, 7:uphi
	diff[0] = var[6];diff[1] = var[4];diff[2] = var[5];diff[3] = var[7];
	double F_r, F_t, F_theta, F_phi;
	int ii, jj;
	double Gamma[4][4][4];
	double u[4]; //Ë³ÐòÊÇut,ur,uth,uphi¡£¡£¡£
	u[0] = var[6];u[1] = var[4];u[2] = var[5];u[3] = var[7];

	Christoffel_KRZ(spin, krz_d, var[1], var[2], Gamma);
	F_r = 0;
	F_theta = 0;
	F_t = 0;
	F_phi = 0;
	for (ii = 0;ii < 4;ii++) {
		for (jj = 0;jj < 4;jj++) {

			F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
			F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
			F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
			F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
		}
	}
	diff[6] = F_t;
	diff[4] = F_r;
	diff[5] = F_theta;
	diff[7] = F_phi;
}
//Ò»¸öÒÔÇ°ÓÃµÄRK4Ëã·¨
//for (i = 0; i < 1000000;i++) {
//	/*if (r > 30) {
//	system("pause");
//	}*/
//	if (var[1] < horizon + 0.01) {
//		printf("Fell into horizon\n");
//		break;
//	}
//	for (;;) {
//		if (var[1] < horizon + 0.01) {
//			printf("Fell into horizon\n");
//			break;
//		}
//		t = var[0];	r = var[1];	th = var[2];	phi = var[3];	ur = var[4];	uth = var[5];	ut = var[6];	uphi = var[7];
//
//		/******************** CALCULATE k1 *****************/
//		metric_KRZ(spin, krz_d, r, th, g);
//		Christoffel_KRZ(spin, krz_d, r, th, Gamma);
//		k1[0] = ut;	varnew1[0] = var[0] + 0.5*dtau*k1[0];
//		k1[3] = uphi;	varnew1[3] = var[3] + 0.5*dtau*k1[3];
//		k1[2] = uth;	varnew1[2] = var[2] + 0.5*dtau*k1[2];
//		k1[1] = ur;		varnew1[1] = var[1] + 0.5*dtau*k1[1];
//
//
//
//		u[0] = ut;
//		u[1] = ur;
//		u[2] = uth;
//		u[3] = uphi;
//
//
//
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		ita = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//				ita += g[ii][jj] * u[ii] * u[jj];
//			}
//		}
//		//if (ita > -0.8) {
//		//system("pause");
//		//break;
//		//ita = -0.8;
//		//}
//		E = -g[0][0] * ut - g[0][3] * uphi;
//		Lz = g[0][3] * ut + g[3][3] * uphi;
//		Q = uth*g[2][2] * uth*g[2][2] + cos(th)*cos(th)*(spin*spin*(mu*mu - E*E) + Lz*Lz / sin(th) / sin(th));
//
//		k1[4] = F_r;	varnew1[4] = var[4] + 0.5*dtau*k1[4];
//		k1[5] = F_theta;	varnew1[5] = var[5] + 0.5*dtau*k1[5];
//		k1[6] = F_t;	varnew1[6] = var[6] + 0.5*dtau*k1[6];
//		k1[7] = F_phi;	varnew1[7] = var[7] + 0.5*dtau*k1[7];
//		/*********************finished calculate k1********************/
//
//		/******************** CALCULATE k2 *****************/
//		metric_KRZ(spin, krz_d, varnew1[1], varnew1[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew1[1], varnew1[2], Gamma);
//
//		ut = varnew1[6];	k2[0] = varnew1[6];	varnew2[0] = var[0] + 0.5*dtau*k2[0];
//		ur = varnew1[4];	k2[1] = varnew1[4];	varnew2[1] = var[1] + 0.5*dtau*k2[1];
//		uth = varnew1[5];	k2[2] = varnew1[5];	varnew2[2] = var[2] + 0.5*dtau*k2[2];
//		uphi = varnew1[7];	k2[3] = varnew1[7];	varnew2[3] = var[3] + 0.5*dtau*k2[3];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k2[4] = F_r;	varnew2[4] = var[4] + 0.5*dtau*k2[4];
//		k2[5] = F_theta;	varnew2[5] = var[5] + 0.5*dtau*k2[5];
//		k2[6] = F_t;	varnew2[6] = var[6] + 0.5*dtau*k2[6];
//		k2[7] = F_phi;	varnew2[7] = var[7] + 0.5*dtau*k2[7];
//
//		/*********************finished calculate k2********************/
//
//		/*********************CALCULATE k3********************/
//		metric_KRZ(spin, krz_d, varnew2[1], varnew2[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew2[1], varnew2[2], Gamma);
//
//		ut = varnew2[6];	k3[0] = varnew2[6];	varnew3[0] = var[0] + dtau*k3[0];
//		ur = varnew2[4];	k3[1] = varnew2[4];	varnew3[1] = var[1] + dtau*k3[1];
//		uth = varnew2[5];	k3[2] = varnew2[5];	varnew3[2] = var[2] + dtau*k3[2];
//		uphi = varnew2[7];	k3[3] = varnew2[7];	varnew3[3] = var[3] + dtau*k3[3];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k3[4] = F_r;	varnew3[4] = var[4] + dtau*k3[4];
//		k3[5] = F_theta;	varnew3[5] = var[5] + dtau*k3[5];
//		k3[6] = F_t;	varnew3[6] = var[6] + dtau*k3[6];
//		k3[7] = F_phi;	varnew3[7] = var[7] + dtau*k3[7];
//
//		/*********************finished calculate k3********************/
//
//
//		/*********************CALCULATE k4********************/
//		metric_KRZ(spin, krz_d, varnew3[1], varnew3[2], g);
//		Christoffel_KRZ(spin, krz_d, varnew3[1], varnew3[2], Gamma);
//
//		ut = varnew3[6];	k4[0] = varnew3[6];
//		ur = varnew3[4];	k4[1] = varnew3[4];
//		uth = varnew3[5];	k4[2] = varnew3[5];
//		uphi = varnew3[7];	k4[3] = varnew3[7];
//
//		u[0] = ut;	u[1] = ur;	u[2] = uth;	u[3] = uphi;
//
//		F_r = 0;
//		F_theta = 0;
//		F_t = 0;
//		F_phi = 0;
//		for (ii = 0;ii < 4;ii++) {
//			for (jj = 0;jj < 4;jj++) {
//
//				F_t -= 0.5*Gamma[0][ii][jj] * u[ii] * u[jj];
//				F_r -= 0.5*Gamma[1][ii][jj] * u[ii] * u[jj];
//				F_theta -= 0.5*Gamma[2][ii][jj] * u[ii] * u[jj];
//				F_phi -= 0.5*Gamma[3][ii][jj] * u[ii] * u[jj];
//			}
//		}
//		k4[4] = F_r;
//		k4[5] = F_theta;
//		k4[6] = F_t;
//		k4[7] = F_phi;
//		/****************************finished calculate k4********************/
//
//		if (fabs(fabs(dtau*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]) / 6.0) - dt) < tol) break;
//		else {
//			dtau = fabs(6.0*dt / (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]));
//		}
//	}
//
