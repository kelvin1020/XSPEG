#ifndef _DEF_H
#include "def.h"
#endif

#include <stdlib.h>
#include <iostream>
#include <string>
//#include <time.h>
using namespace std;
#ifndef max
#define max(a,b) (( (a) >= (b)) ? (a) : (b))
#endif
#define N 16 //8个变量
#define TIMELENGTH 1000000
double orbitlist[N][TIMELENGTH] = { 0 };
//double evolv_e[TIMELENGTH] = { 0 }, evolv_p[TIMELENGTH] = { 0 }, evolv_iota[TIMELENGTH] = { 0 };
double current_e, current_p, current_iota, current_omgr, current_omgth, current_omgphi;
int lengthflag = 0;//0表示算出来的轨道还太短，先用输入的e, p, iota
int timeindex = 0;
double massratio;

void equations(double var[],double diff[],double spin,double defpar[],double massratio);
void dboydcar(double spin, double car[], double boy[], double dboy_dcar[][4]);
void dboydcar_boyknown(double spin, double car[], double boy[], double dboy_dcar[][4]);
void dcardboy(double spin, double boy[], double car[], double dcar_dboy[][4]);
void highderiv(double spin, double defpar[], double boy[], double boy_taudot[], double x[][4]);
void radacc(double spin, double x[][4], double acc[]);
void Gamma_car(double spin, double defpar[], double car[], double Gamma[][4][4]);
void get_current_orbitpar(double tlist[], double rlist[], double thlist[], double philist[], double orbitpar[]);
int time_long_enough(double rlist[], double thlist[]);
void get_flux(double spin,double E, double L, double Q, double flux[]);

int main(int argc, char *argv[])
{//RK45的参数
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

	double spin, spin2, defpar[10] = { 0 };
	double isco, xin;
	double robs_i, robs_f;
	double pstep;
	int i, j, k, m;
	int ii,jj;
	char filename_o[128];

	FILE *foutput, *finput;

	double r, th, t, phi, ur, ut, uth, uphi, tau=0, dtau; //坐标，4速度，固有时，固有时步长

	//double rnew, thnew, tnew, phinew, urnew, utnew, uthnew, uphinew, dtau;

	double E, Lz;

	double F_theta, F_r,F_t,F_phi;

	double g[4][4] = { 0 };
	double Gamma[4][4][4] = { 0 };
	double u[4] = { 0 };
	double k1[N], k2[N], k3[N], k4[N],k5[N],k6[N];//16个变量的顺序是t,r,\theta,\phi,u^t, u^r,u^\theta, u^\phi, ts,rs,\theta s,\phi s,us^t, us^r,us^\theta, us^\phi
	double var[N];
	double ita,mu=-1;//计算过程中四速度的模是\ita， 描述粒子性质的是mu 有质量：mu=-1，无质量： mu=0

//input argument: 1. M, 2. spin, 3. E, 4. Lz, 5. Q, 6. r0, 7. tottime, 8-10. defpar, 11. mass ratio
	
	/***************************改成秒间隔需要改的部分********************/

	double dt;
	double Grav, clight, Msol, M, dt_sec = 0.1;
	double t_sec = 0;
	
	Grav = 6.674e-11;//#引力常数
	clight = 2.998e8;//#光速
	Msol = 1.989e30;  //#太阳质量，以千克做单位
	
	M = atof(argv[1]);
	//dt_sec = dt*M*Msol*Grav / clight / clight / clight

	dt = dt_sec*clight*clight*clight / M / Msol / Grav;
	
	/***************************改成秒间隔需要改的部分********************/

	dtau = dt;//初始步长
	double tottime=atof(argv[7]);
	massratio = atof(argv[11]);//atof(argv[11]);

	current_e = atof(argv[12]);
	current_p = atof(argv[13]);
	current_iota = atof(argv[14]);

	spin = atof(argv[2]);//+0.001*0.5;
	//for (spin = 0.55428452790227312;spin <  0.55428452790227312 +0.01;spin += 0.05) {
		printf("Metric info: spin=%f M=%f\n", spin,M);
		defpar[1] = atof(argv[8]),defpar[2]=atof(argv[9]), defpar[3]=atof(argv[10]);//+0.001*0.2;
//		for (defpar[2] = 0.0; defpar[2] <= 0.0; defpar[2] += 0.1) {
			printf("deformation parameters: d1=%f, d2=%f, d3=%f\n",defpar[1], defpar[2],defpar[3]);
			t = 0;
			//r = 13.0;
			th = Piby2;
			phi = 0;


			//ut = ( E*g[3][3]+Lz*g[0][3] ) / ( g[0][3]*g[0][3] - g[0][0]*g[3][3] );
			//uphi = (E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]) ;
			//uth = sqrt((-1 - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));

			/********************↓圆轨道的能量和角动量 from https://arxiv.org/pdf/1105.2959.pdf （好像是错的。。）↓******************/
			/*	Lz = fabs(spin*spin + 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) - 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, corotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin - sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) - 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, corotating

			Lz = fabs(spin*spin - 2 * spin*sqrt(r) + r*r) / sqrt(r*r*(r - 3) + 2 * spin*sqrt(r*r*r));//Kerr_circular_orbit, counterrotating
			E = ((pow(r, 1.25)*fabs((spin*spin + r*(r - 2))*(spin + sqrt(r*r*r))) / sqrt((r - 3)*sqrt(r) + 2 * spin)) + 2 * spin*r*Lz) / (r*(r*r*r + spin*spin*(r + 2))); //Kerr_circular_orbit, counterrotating
			*/
			/********************↑from https://arxiv.org/pdf/1105.2959.pdf （好像是错的。。。）↑******************/

			/********************↓圆轨道的能量和角动量 from http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf （这个是对的）↓******************/
			//Lz = ( sqrt(r) -2*spin/r + spin*spin/sqrt(r*r*r) ) / sqrt(1 - 3 / r + 2 * spin / sqrt(r*r*r));
			//E = (1-2/r + spin /sqrt(r*r*r) ) / sqrt(1-3/r + 2*spin/sqrt(r*r*r) );
			/********************↑from http://www.tapir.caltech.edu/~chirata/ph236/2011-12/lec27.pdf ↑******************/


			double horizon, r0, th0 = Piby2, phi0 = 0, t0 = 0;//初始条件
			double ecc, p;//eccentricity和rectum 
			double rmax, rmin, invgmax[4][4], invgmin[4][4];//r的上下限以及该处的g
			double invg[4][4];
			double EoverL;//由e和p算E和Lz的中间变量
			double EoverL2, E2, L2;
			//double iota = Pi / 6;
			double Q,Q0;//carter constant
			horizon = 1 + sqrt(1 - spin*spin);
			/***************************赤道面上由e,p决定的轨道↓******************/

			//for (ecc = 0.41149551640029858;ecc <= 0.41149551640029858;ecc = ecc + 0.004) {
			/*for (ecc = 0.5;ecc <= 0.5;ecc = ecc + 0.004) {
				//	for (p = 6.4825501282607396;p <= 6.4825501282607396;p = p + 0.04) {
				for (p = 6;p <= 6;p = p + 0.04) {*/
					
			//ecc = 0.43173473149300562,p= 7.1717260434110983;
					E = atof(argv[3]), Lz=atof(argv[4]), Q=atof(argv[5]);//python解出来的   and note that it's not useful actually, Q is calculated below, and for non-Kerr case this is just initial Q
					th0 = Piby2;//赤道面上carter constant和theta方向速度关系比较简单
					r0 = atof(argv[6]);//看看取初始在rmax会怎么样？
					Q0 = Q;
			/*
					rmax = p / (1 - ecc);
					rmin = p / (1 + ecc);
					r0 = rmax;
					metric_inverse(spin, defpar, rmax, th0, invgmax);
					metric_inverse(spin, defpar, rmin, th0, invgmin);


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
					metric(spin, defpar, r, th, g);
					metric_inverse(spin, defpar, r, th, invg);



					/*
					spin = 0.900000;current_p = 6.000000;current_e = 0.500000;current_iota = 0.700671;
					E = 0.9425685155;Lz = 2.2583745750, Q = 3.6658018921;
					double myflux[3];
					massratio = 1;
					get_flux(spin, Lz, Q, myflux);
					printf("%.10e %.10e %.10e", myflux[0], myflux[1], myflux[2]);
					system("pause");
					*/
					printf("Mass ratio is %f\n",massratio);
					printf("Starting e=%.6f p=%.6f iota=%.6f\n\n",  current_e, current_p, current_iota);
					/********************↓二分找E,Lz↓******************/
					/*
					if (abs(ecc - 0.0) < 1e-6) {
						r = p;r0 = p;
						double upE = 1, downE = 0.9, eps = 1e-10;
						double invg[4][4];
						double curf, upf, downf;
						Christoffel(spin, defpar, r, th, Gamma);
						metric_inverse(spin, defpar, r, th, invg);

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

					/********************↑二分找圆轨道 ↑******************/
					/*
					double testE, testL, testut, testup, utoverup;
					if (abs(ecc - 0.0) < 1e-6) {
						///解方程找圆轨道↓
						//p = 10;
						Christoffel(spin, defpar, p, th, Gamma);
						metric(spin, defpar, p, th, g);
						utoverup = (-Gamma[1][0][3] + sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
						testup = sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
						testut = testup*utoverup;
						testE = -g[0][0] * testut - g[0][3] * testup;
						testL = g[0][3] * testut + g[3][3] * testup;
						E = testE;
						Lz = testL;
						//   解方程找圆轨道↑
					}
					*/

					/*looking for circular orbit↓*/
					/*
					metric(spin, defpar, r, th, g);
					Christoffel(spin, defpar, r, th, Gamma);
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

					/*looking for circular orbit↑*/

					ut = (E*g[3][3] + Lz*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					uphi = -(E*g[0][3] + Lz*g[0][0]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
					//uth = sqrt(Q) / g[2][2];
					//ur = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[2][2] * uth*uth) / (g[1][1]));
					ur = 0;
					if (abs((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2])) < eps) uth = 0;
					else uth = sqrt((mu - g[0][0] * ut*ut - 2 * g[0][3] * ut*uphi - g[3][3] * uphi*uphi - g[1][1] * ur*ur) / (g[2][2]));
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

					//sprintf(filename_o, "circular_trace_spin%.2f_d%.2f_r%.2f.dat", spin, defpar,r);


					//测试坐标变换
					/*
					double car[4]; double boy[4];
					double dboy_dcar[4][4]; double dcar_dboy[4][4];

					boy[0] = t;boy[1] = r;boy[2] = th;boy[3] = phi;
					dcardboy(spin, boy, car, dcar_dboy);
					cout << car[0]<<car[1]<<car[2]<<car[3];

					dboydcar(spin, car, boy, dboy_dcar);
					int flag = 0;
					for (int ind = 0;ind < 4;ind++) {
						if (!isfinite(boy[ind])) {
							flag = 1;
							cout << "boy [" << ind << ']' << endl;
						}
						if (!isfinite(car[ind])) {
							flag = 1;
							cout << "car [" << ind << ']' << endl;
						}
						for (int ind2 = 0;ind2 < 4;ind2++) {
							if (!isfinite(dboy_dcar[ind][ind2])) {
								flag = 1;
								cout << "dboy_dcar [" << ind << "]["<<ind2<<']' << endl;
							}
							if (!isfinite(dcar_dboy[ind][ind2])) {
								flag = 1;
								cout << "dcar_dboy [" << ind << "]["<<ind2<<']' << endl;
							}
						}
					}
					if (flag == 1) system("pause");*/
					//上面在测试坐标变换

					sprintf(filename_o,"ORBCAR");// "trace_M%.0f_spin%.6f_E%.6f_Lz%.6f_Q%.6f_d1%.6f_d2%.6f_d3%.6f.dat", M, spin, E, Lz, Q, defpar[1], defpar[2], defpar[3]);
					foutput = fopen(filename_o, "w");
					/*
					FILE *foutputnk;
					sprintf(filename_o, "NK_M%.0f_spin%.6f_E%.6f_Lz%.6f_Q%.6f_d1%.6f_d2%.6f_d3%.6f.dat", M, spin, E, Lz, Q, defpar[1], defpar[2], defpar[3]);
					foutputnk = fopen(filename_o, "w");
					FILE *foutputflux;
					sprintf(filename_o, "flux_M%.0f_spin%.6f_E%.6f_Lz%.6f_Q%.6f_d1%.6f_d2%.6f_d3%.6f.dat", M, spin, E, Lz, Q, defpar[1], defpar[2], defpar[3]);
					foutputflux = fopen(filename_o, "w");
					*/
					//Zprintf("index\t tau\t t\t r\t theta\t phi\t ut\t ur\t uth\t uphi\t ita\t E\t Lz\n");

					var[0] = t;	var[1] = r;	var[2] = th;	var[3] = phi;	var[4] = ut;	var[5] = ur;	var[6] = uth;	var[7] = uphi;
					var[8] = 0;	var[9] = 0;	var[10] = 0;	var[11] = 0;	var[12] = 0;	var[13] = 0;	var[14] = 0;	var[15] = 0;//小量初始为0
					for (i = 0; i < N; i++) {
						orbitlist[timeindex][i] = var[i];
					}

					double E0 = E, L0 = Lz;
					double h;//RK45, adaptive step
					double diff[N];
					double vars_temp[N];//temp used in RK45
					double z[N], y[N];//used to estimate error
					double err, maxerr=1e-10, minerr=1e-11;//计算过程中的误差，容许的最大误差，容许的最小误差
					int check=0;//标记有没有超过最大误差或小于最小误差
					int index = 0;
					h = dtau;
					tau = 0;
					double * tlist, *rlist, *thlist, *philist;
					double orbitpar[6];
					double percentage=0;//time_t optime = time(NULL);
					for(;t_sec<tottime;){

						if (var[1] + var[9] < horizon*1.1) {
							printf("Plunge at r=%f\n", var[1] + var[9]);
							
                                                        //printf("\033[1A"); // move cursor one line up
					                //printf("\033[1A"); // move cursor one line up
        		                                //printf("\033[K");   // delete till end of line
	
							//system("pause");
							break;
						}
						//if(t_sec/tottime*100>percentage+1.0){
                                                        percentage=t_sec/tottime*100;
                                                        printf("\033[1A"); // move cursor one line up
                                                        printf("\033[K");   // delete till end of line

                                                        printf("Time t/M=%f. Current semilatus p=%f eccentricity e=%f inclination iota=%f \n", var[0]+var[8],current_p,current_e,current_iota);
                                                //}

//						printf("p=%f\n", current_p);
						for (i = 0;i < 8;i++) {
							if (abs(var[i + 8]) > 1e-4*abs(var[i])) {
								var[i] = var[i] + var[i + 8];
								var[i + 8] = 0;
							}
						}

						tlist = orbitlist[0];
						rlist = orbitlist[1];
						thlist = orbitlist[2];
						philist = orbitlist[3];
						
						if (lengthflag == 0) {
							if (time_long_enough(rlist, thlist)) {
								lengthflag = 1;
							}
							//if (timeindex>5&&current_e < 1e-4) lengthflag = 1;
						}
						
						
						t_sec = var[0] * M*Msol*Grav / clight / clight / clight;

						metric(spin, defpar, var[1]+var[9], var[2]+var[10], g);
						ut = var[4]+var[12];ur = var[5]+var[13];uth = var[6]+var[14];uphi = var[7]+var[15];
						
						ita = g[0][0] * ut*ut + g[1][1] * ur*ur + g[2][2] * uth*uth + g[3][3] * uphi*uphi + 2 * g[3][0] * uphi*ut;
						E = -g[0][0] * ut - g[0][3] * uphi;
						Lz = g[0][3] * ut + g[3][3] * uphi;
						Q = uth*g[2][2] * uth*g[2][2] + cos(var[2]+var[10])*cos(var[2] + var[10])*(spin*spin*(ita*ita - E*E) + Lz*Lz / sin(var[2] + var[10]) / sin(var[2] + var[10]));
						
						//下面在测试坐标变换
						/*
						double car[4]; double boy[4];
						double dboy_dcar[4][4]; double dcar_dboy[4][4];
						boy[0] = var[0];boy[1] = var[1];boy[2] = var[2];boy[3] = car[3];
						dcardboy(spin, boy, car, dcar_dboy);
						cout << "car:   "<< car[0] <<'\t'<< car[1] <<'\t'<< car[2] <<'\t'<< car[3]<<endl;

						dboydcar(spin, car, boy, dboy_dcar);
						int flag = 0;
						for (int ind = 0;ind < 4;ind++) {
							if (!isfinite(boy[ind])) {
								flag = 1;
								cout << "boy [" << ind << ']' << endl;
							}
							if (!isfinite(car[ind])) {
								flag = 1;
								cout << "car [" << ind << ']' << endl;
							}
							for (int ind2 = 0;ind2 < 4;ind2++) {
								if (!isfinite(dboy_dcar[ind][ind2])) {
									flag = 1;
									cout << "dboy_dcar [" << ind << "][" << ind2 << ']' << endl;
								}
								if (!isfinite(dcar_dboy[ind][ind2])) {
									flag = 1;
									cout << "dcar_dboy [" << ind << "][" << ind2 << ']' << endl;
								}
							}
						}
						if (flag == 1) system("pause");*/
						//上面在测试坐标变换

						check = 0;
						equations(var, diff,spin,defpar,massratio);
						for (i = 0; i < N; i++)
						{
							k1[i] = h*diff[i];
							vars_temp[i] = var[i] + a1*k1[i];
						}

						equations(vars_temp, diff,spin,defpar,massratio);
						for (i = 0; i < N; i++)
						{
							k2[i] = h*diff[i];
							vars_temp[i] = var[i] + b1*k1[i] + b2*k2[i];
						}

						equations(vars_temp, diff,spin,defpar, massratio);
						for (i = 0; i < N; i++)
						{
							k3[i] = h*diff[i];
							vars_temp[i] = var[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
						}

						equations(vars_temp, diff,spin,defpar, massratio);
						for (i = 0; i < N; i++)
						{
							k4[i] = h*diff[i];
							vars_temp[i] = var[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
						}

						equations(vars_temp, diff,spin,defpar, massratio);
						for (i = 0; i < N; i++)
						{
							k5[i] = h*diff[i];
							vars_temp[i] = var[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
						}

						equations(vars_temp, diff,spin,defpar, massratio);
						for (i = 0; i < N; i++)
							k6[i] = h*diff[i];

						for (i = 0; i < N; i++)
						{
							y[i] = var[i] + x1*k1[i] + x2*k2[i] + x3*k3[i] + x4*k4[i] + x5*k5[i];
							z[i] = var[i] + z1*k1[i] + z2*k2[i] + z3*k3[i] + z4*k4[i] + z5*k5[i] + z6*k6[i];
							if (i<N/2 && abs(var[i] - 0) < eps) continue;
							err = fabs((y[i] - z[i]) / max(fabs(var[i]), fabs(y[i]) ) );
							if (err > maxerr) {
								check = 1;//有些误差太大
							}
							else if (err < minerr&&check != 1) {
								check = -1;//所有误差都太小
							}
						}
						if (h < 0.1) {
							check = -1;
						}
						
						//if (index == 565) system("pause");
						for (i = 0;i < N;i++) {
							if (!isfinite(y[i])) {
								printf("Numerical Error\n");
								break;//system("pause");
							}
						}
						if (check == 1) {
							h /= 1.1;
						}
						
						else if (check == -1) {
							for (i = 0;i < N / 2;i++) {
								k1[i] = k1[i] + k1[i + 8];
								k2[i] = k2[i] + k2[i + 8];
								k3[i] = k3[i] + k3[i + 8];
								k4[i] = k4[i] + k4[i + 8];
								k5[i] = k5[i] + k5[i + 8];
								k6[i] = k6[i] + k6[i + 8];
							}
							
							fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n",
								index, t_sec, tau, var[0]+var[8], var[1] + var[9], var[2] + var[10], var[3] + var[11], var[4] + var[12], var[5] + var[13], var[6] + var[14], var[7] + var[15],
								/*F_t*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]),
								/*F_r*/(z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5]),
								/*F_theta*/z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6],
								/*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
							//printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n",
							//	index, t_sec, t, var[0], var[1] + var[9], var[2] + var[10], var[3] + var[11], var[4] + var[12]/*ut*/, var[5] + var[13]/*ur*/, var[6] + var[14]/*uth*/, var[7] + var[15]/*uphi*/,
							//	mu, E, Lz, Q);
							//fprintf(foutputnk, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
							//	var[0]+var[8], var[1] + var[9], var[2] + var[10], var[3] + var[11], E, Lz, Q);
							fflush(foutput);
							//fflush(foutputnk);
							//double flux[3];
							//get_flux(spin, Lz, Q, flux);
							//printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", var[0] + var[8], var[1]+var[9], current_e, current_p, current_iota, E, Lz, Q, ita, flux[0], flux[1], flux[2]);
							/*
							fprintf(foutputflux, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
								tau, var[0] + var[8],ut, E, Lz, Q,flux[0],flux[1],flux[2]);
							fflush(foutputflux);*/
							tau = tau + h;
							h *= 1.1;
							index++;
							timeindex++;
							for (i = 0; i < N; i++) {
								var[i] = y[i];
								if (i<N/2) orbitlist[i][timeindex] = var[i]+var[i+N/2];
							}
							if (lengthflag == 1) {//更新轨道参数

								get_current_orbitpar(tlist, rlist, thlist, philist, orbitpar);
								current_e = orbitpar[0];
								current_p = orbitpar[1];
								current_iota = orbitpar[2];
								current_omgr = orbitpar[3];
								current_omgth = orbitpar[4];
								current_omgphi = orbitpar[5];
								//printf("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", orbitpar[0], orbitpar[1], orbitpar[2], orbitpar[3], orbitpar[4], orbitpar[5]);

							}
						}
						else {
							for (i = 0;i < N / 2;i++) {
								k1[i] = k1[i] + k1[i + 8];
								k2[i] = k2[i] + k2[i + 8];
								k3[i] = k3[i] + k3[i + 8];
								k4[i] = k4[i] + k4[i + 8];
								k5[i] = k5[i] + k5[i + 8];
								k6[i] = k6[i] + k6[i + 8];
							}
							fprintf(foutput, "%d\t %.10f\t %.10f\t %.10f\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e\t %.10e \n",
								index, t_sec, tau, var[0] + var[8], var[1] + var[9], var[2] + var[10], var[3] + var[11], var[4] + var[12], var[5] + var[13], var[6] + var[14], var[7] + var[15],
								/*F_t*/(z1*k1[4] + z2*k2[4] + z3*k3[4] + z4*k4[4] + z5*k5[4] + z6*k6[4]),
								/*F_r*/(z1*k1[5] + z2*k2[5] + z3*k3[5] + z4*k4[5] + z5*k5[5] + z6*k6[5]),
								/*F_theta*/z1*k1[6] + z2*k2[6] + z3*k3[6] + z4*k4[6] + z5*k5[6] + z6*k6[6],
								/*F_phi*/z1*k1[7] + z2*k2[7] + z3*k3[7] + z4*k4[7] + z5*k5[7] + z6*k6[7]);
							//printf("%d\t %.10f\t %.3f\t %.3f\t %.6f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.10e \t %.10e \t %.10e \t %.10e\n",
							//	index, t_sec, t, var[0], var[1] + var[9], var[2] + var[10], var[3] + var[11], var[4] + var[12]/*ut*/, var[5] + var[13]/*ur*/, var[6] + var[14]/*uth*/, var[7] + var[15]/*uphi*/,
							//	mu, (E - E0) / E0, (Lz - L0) / L0, Q);
							//fprintf(foutputnk, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
							//	var[0] + var[8], var[1] + var[9], var[2] + var[10], var[3] + var[11], E, Lz, Q);
							fflush(foutput);
							//fflush(foutputnk);
							//double flux[3];
							//get_flux(spin, Lz, Q, flux);
							//printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", var[0] + var[8],var[1]+var[9],current_e,current_p,current_iota, E, Lz, Q, ita, flux[0], flux[1], flux[2]);
							/*
							fprintf(foutputflux, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
								tau, var[0] + var[8], ut, E, Lz, Q, flux[0], flux[1], flux[2]);
							fflush(foutputflux);*/
							tau = tau + h;
							index++;
							timeindex++;
							for (i = 0; i < N; i++) {
								var[i] = y[i];
								if (i<N / 2) orbitlist[i][timeindex] = var[i] + var[i + N / 2];

							}
							if (lengthflag == 1) {//更新轨道参数

								get_current_orbitpar(tlist, rlist, thlist, philist, orbitpar);
								current_e = orbitpar[0];
								current_p = orbitpar[1];
								current_iota = orbitpar[2];
								current_omgr = orbitpar[3];
								current_omgth = orbitpar[4];
								current_omgphi = orbitpar[5];
								//printf("%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", orbitpar[0], orbitpar[1], orbitpar[2], orbitpar[3], orbitpar[4], orbitpar[5]);

							}
						}

					}
					fclose(foutput);
					
                                        //printf("\033[1A"); // move cursor one line up
                                        //printf("\033[K");   // delete till end of line
                                        if (var[1] + var[9] > horizon*1.1) printf("reach time limit. \n");


//time_t edtime = time(NULL);
					//cout << (edtime - optime) << endl;
					//system("pause");
				
			
		//}

	//}
	return 0;

}

void equations(double var[], double diff[],double spin,double defpar[],double massratio) {//由8个变量算斜率
	diff[0] = var[4];diff[1] = var[5];diff[2] = var[6];diff[3] = var[7];
	
	int ii, jj;
	double Gamma[4][4][4];
	double u[4]; 
	u[0] = var[4]+var[12];u[1] = var[5]+var[13];u[2] = var[6]+var[14];u[3] = var[7]+var[15];

	Christoffel(spin, defpar, var[1]+var[9], var[2]+var[10], Gamma);
	diff[4] = 0;
	diff[5] = 0;
	diff[6] = 0;
	diff[7] = 0;
	for (ii = 0;ii < 4;ii++) {
		for (jj = 0;jj < 4;jj++) {

			diff[4] -= Gamma[0][ii][jj] * u[ii] * u[jj];
			diff[5] -= Gamma[1][ii][jj] * u[ii] * u[jj];
			diff[6] -= Gamma[2][ii][jj] * u[ii] * u[jj];
			diff[7] -= Gamma[3][ii][jj] * u[ii] * u[jj];
		}
	}
	
	//检验highderiv
	/*
	double x[8][4] = { 0 }, boy[4];
	boy[0] = var[0]+var[8];boy[1] = var[1]+var[9]; boy[2] = var[2]+var[10]; boy[3] = var[3]+var[11];
	highderiv(spin, defpar, boy, u, x);
	int flag = 0;
	for (int ind = 0;ind < 8;ind++) {
		
		for (int ind2 = 0;ind2 < 4;ind2++) {
			if (!isfinite(x[ind][ind2])) {
				flag = 1;
				cout << "x [" << ind << "][" << ind2 << ']' << endl;
			}
			
		}
	}
	if (flag == 1) system("pause");
	*/
	//上面在检验highderiv

	diff[8] = var[12];
	diff[9] = var[13];
	diff[10] = var[14];
	diff[11] = var[15];
	/*
	double x[8][4] = { 0 }, boy[4], acc[4] = { 0 }, acc_car[4] = { 0 };
	boy[0] = var[0] + var[8];boy[1] = var[1] + var[9]; boy[2] = var[2] + var[10]; boy[3] = var[3] + var[11];
	highderiv(spin, defpar, boy, u, x);
	radacc(spin, x, acc_car);
	for (int i = 1;i < 4;i++) acc_car[i] = acc_car[i] * massratio*u[0] * u[0] - var[i + 12] * F_t/u[0];//乘上质量比才是真*加速度，另外要转换成对tau的加速度
	double car[4] = { 0 }, dboy_dcar[4][4] = { 0 };
	car[0] = boy[0];
	car[1] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * cos(boy[3]);
	car[2] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * sin(boy[3]);
	car[3] = boy[1] * cos(boy[2]);
	dboydcar_boyknown(spin, car, boy, dboy_dcar);
	for (int mu = 0;mu < 4;mu++) {
		for (int al = 0;al < 4;al++) {
			acc[mu]+= acc_car[al] * dboy_dcar[mu][al];
		}
	}*/
	double F_r=0, F_t=0, F_theta=0, F_phi=0;//self-force, 注意其实是上标
	double flux[3] = { 0 };
	double r = var[1] + var[9], th = var[2] + var[10];
	double E, Lz, Q, ut = u[0], ur = u[1], uth = u[2], uphi = u[3];
	double g[4][4] = { 0 };
	metric(spin, defpar, r, th, g);
	E = -g[0][0] * ut - g[0][3] * uphi;
	Lz = g[0][3] * ut + g[3][3] * uphi;
	Q = uth*g[2][2] * uth*g[2][2] + cos(th)*cos(th)*(spin*spin*(1 - E*E) + Lz*Lz / sin(th) / sin(th));
	get_flux(spin,E, Lz, Q, flux);
	if (r < 7) {
		r=r;
	}

	double Edot = flux[0] * ut, Ldot = flux[1] * ut, Qdot = flux[2] * ut;//对tau的flux

	//Qdot = 0;

	if (abs(g[0][3] * g[0][3] - g[0][0] * g[3][3]) < 1e-8) {
		F_t = 0;F_phi = 0;
	}
	else {
		F_t = (Edot*g[3][3] + Ldot*g[0][3]) / (g[0][3] * g[0][3] - g[0][0] * g[3][3]);
		F_phi = (Edot*g[0][3] + Ldot*g[0][0]) / ( g[0][0] * g[3][3] - g[0][3] * g[0][3]);
	}
	if (abs(uth) < 1e-3) F_theta = 0;
	else {
		F_theta = ( Qdot - (2 * cos(th) * cos(th) *spin*spin *E*Edot + 2 * cos(th)*cos(th) / sin(th) / sin(th) * Lz *Ldot) )/2.0/g[2][2]/g[2][2]/uth;
	}

	//F_theta = 0;

	if (abs(ur) < 1e-1*massratio || abs(ur)<1e-16 ) F_r = 0;
	else {
		F_r = -(g[0][0] * ut*F_t + g[0][3] * ut*F_phi + g[0][3] * uphi*F_t + g[2][2] * uth*F_theta + g[3][3] * uphi*F_phi) / g[1][1] / ur;
	}
	//F_r = 0;
	diff[12] = F_t;
	diff[13] = F_r;
	diff[14] = F_theta;
	diff[15] = F_phi;
}
int time_long_enough(double rlist[], double thlist[]) {
	double ra[2], rp[2], thmin[2], thmax[2];
	int tind_ra[2], tind_rp[2], tind_thmin[2], tind_thmax[2];
	int rai = 0, rpi = 0, thmini = 0, thmaxi = 0;
	int tind;
	if (current_e > 1e-4&&current_iota>1e-4) {
		for (tind = timeindex - 1;tind >= 1;tind--) {
			if (rlist[tind] > current_p && rlist[tind] >= rlist[tind + 1] && rlist[tind] >= rlist[tind - 1] && rai < 2) {
				ra[rai] = rlist[tind];
				tind_ra[rai] = tind;
				rai++;
			}
			if (rlist[tind] < current_p && rlist[tind] <= rlist[tind + 1] && rlist[tind] <= rlist[tind - 1] && rpi < 2) {
				rp[rpi] = rlist[tind];
				tind_rp[rpi] = tind;
				rpi++;
			}
			if (thlist[tind] > Piby2 &&thlist[tind] >= thlist[tind + 1] && thlist[tind] >= thlist[tind - 1] && thmaxi < 2) {
				thmax[thmaxi] = thlist[tind];
				tind_thmax[thmaxi] = tind;
				thmaxi++;
			}
			if (thlist[tind] < Piby2 && thlist[tind] <= thlist[tind + 1] && thlist[tind] <= thlist[tind - 1] && thmini < 2) {
				thmin[thmini] = thlist[tind];
				tind_thmin[thmini] = tind;
				thmini++;
			}
			if (rai >= 2 && rpi >= 2 && thmini >= 2 && thmaxi >= 2) break;

		}
	}
	else if (current_iota>1e-4) {
		ra[0] = rlist[timeindex];rai = 1;
		rp[0] = rlist[timeindex];rpi = 1;
		for (tind = timeindex - 1;tind >= 1;tind--) {

			if (thlist[tind] > Piby2 &&thlist[tind] >= thlist[tind + 1] && thlist[tind] >= thlist[tind - 1] && thmaxi < 2) {
				thmax[thmaxi] = thlist[tind];
				tind_thmax[thmaxi] = tind;
				thmaxi++;
			}
			if (thlist[tind] < Piby2 && thlist[tind] <= thlist[tind + 1] && thlist[tind] <= thlist[tind - 1] && thmini < 2) {
				thmin[thmini] = thlist[tind];
				tind_thmin[thmini] = tind;
				thmini++;
			}
			if (thmini >= 2 && thmaxi >= 2) break;

		}
	}
	else if (current_e > 1e-4) {
		thmin[0] = thlist[timeindex];thmini = 1;
		thmax[0] = thlist[timeindex];thmaxi = 1;
		for (tind = timeindex - 1;tind >= 1;tind--) {
			if (rlist[tind] > current_p && rlist[tind] >= rlist[tind + 1] && rlist[tind] >= rlist[tind - 1] && rai < 2) {
				ra[rai] = rlist[tind];
				tind_ra[rai] = tind;
				rai++;
			}
			if (rlist[tind] < current_p && rlist[tind] <= rlist[tind + 1] && rlist[tind] <= rlist[tind - 1] && rpi < 2) {
				rp[rpi] = rlist[tind];
				tind_rp[rpi] = tind;
				rpi++;
			}
			if (rai >= 2 && rpi >= 2) break;

		}
	}
	else {
		ra[0] = rlist[timeindex];rai = 1;
		rp[0] = rlist[timeindex];rpi = 1;
		thmin[0] = thlist[timeindex];thmini = 1;
		thmax[0] = thlist[timeindex];thmaxi = 1;
	}
	if (rai >= 1 && rpi >= 1 && thmini >= 1 && thmaxi >= 1) {
		return 1;
	}
	else return 0;


}

void get_current_orbitpar(double tlist[], double rlist[], double thlist[], double philist[], double orbitpar[]) {
	if (lengthflag == 0) return;

	double ra[2], rp[2], thmin[2], thmax[2];
	int tind_ra[2], tind_rp[2], tind_thmin[2], tind_thmax[2];
	int rai = 0, rpi = 0, thmini = 0, thmaxi = 0;
	int tind;
	if (current_e > 1e-4&&current_iota>1e-4) {
		for (tind = timeindex - 1;tind >= 1;tind--) {
			if (rlist[tind] > current_p && rlist[tind] >= rlist[tind + 1] && rlist[tind] >= rlist[tind - 1] && rai < 2) {
				ra[rai] = rlist[tind];
				tind_ra[rai] = tind;
				rai++;
			}
			if (rlist[tind] < current_p && rlist[tind] <= rlist[tind + 1] && rlist[tind] <= rlist[tind - 1] && rpi < 2) {
				rp[rpi] = rlist[tind];
				tind_rp[rpi] = tind;
				rpi++;
			}
			if (thlist[tind] > Piby2 &&thlist[tind] >= thlist[tind + 1] && thlist[tind] >= thlist[tind - 1] && thmaxi < 2) {
				thmax[thmaxi] = thlist[tind];
				tind_thmax[thmaxi] = tind;
				thmaxi++;
			}
			if (thlist[tind] < Piby2 && thlist[tind] <= thlist[tind + 1] && thlist[tind] <= thlist[tind - 1] && thmini < 2) {
				thmin[thmini] = thlist[tind];
				tind_thmin[thmini] = tind;
				thmini++;
			}
			if (rai >= 2 && rpi >= 2 && thmini >= 2 && thmaxi >= 2) break;

		}
	}
	else if(current_iota>1e-4){
		ra[0] = rlist[timeindex];rai = 1;
		rp[0] = rlist[timeindex];rpi = 1;
		for (tind = timeindex - 1;tind >= 1;tind--) {
			
			if (thlist[tind] > Piby2 &&thlist[tind] >= thlist[tind + 1] && thlist[tind] >= thlist[tind - 1] && thmaxi < 2) {
				thmax[thmaxi] = thlist[tind];
				tind_thmax[thmaxi] = tind;
				thmaxi++;
			}
			if (thlist[tind] < Piby2 && thlist[tind] <= thlist[tind + 1] && thlist[tind] <= thlist[tind - 1] && thmini < 2) {
				thmin[thmini] = thlist[tind];
				tind_thmin[thmini] = tind;
				thmini++;
			}
			if (thmini >= 2 && thmaxi >= 2) break;

		}
	}
	else if (current_e > 1e-4) {
		thmin[0] = thlist[timeindex];thmini = 1;
		thmax[0] = thlist[timeindex];thmaxi = 1;
		for (tind = timeindex - 1;tind >= 1;tind--) {
			if (rlist[tind] > current_p && rlist[tind] >= rlist[tind + 1] && rlist[tind] >= rlist[tind - 1] && rai < 2) {
				ra[rai] = rlist[tind];
				tind_ra[rai] = tind;
				rai++;
			}
			if (rlist[tind] < current_p && rlist[tind] <= rlist[tind + 1] && rlist[tind] <= rlist[tind - 1] && rpi < 2) {
				rp[rpi] = rlist[tind];
				tind_rp[rpi] = tind;
				rpi++;
			}
			if (rai >= 2 && rpi >= 2) break;

		}
	}
	else {

		ra[0] = rlist[timeindex];rai = 1;
		rp[0] = rlist[timeindex];rpi = 1;
		thmin[0] = thlist[timeindex];thmini = 1;
		thmax[0] = thlist[timeindex];thmaxi = 1;
	}
	double realra, realrp, realthmin, realthmax;
	if (rai >= 2) {
		realra = max((ra[0]), (ra[1]));
	}
	else realra = ra[0];
	if (rpi >= 2) {
		realrp = max((rp[0]), (rp[1]));
	}
	else realrp = rp[0];
	if (thmini >= 2) {
		realthmin = max((thmin[0]), (thmin[1]));
	}
	else realthmin = thmin[0];
	if (thmaxi >= 2) {
		realthmax = max((thmax[0]), (thmax[1]));
	}
	else realthmax = thmax[0];

	orbitpar[0] = (realra - realrp) / (realra + realrp);//e
	orbitpar[1] = 2 * realra*realrp / (realra + realrp);//p
	orbitpar[2] = max((realthmax - Piby2), (Piby2 - realthmin));//iota
	/*
	orbitpar[3] = 0;orbitpar[4] = 0;orbitpar[5] = 0;
	if (rai >= 2) {
		orbitpar[3] = 2 * Pi / (tlist[tind_ra[0]] - tlist[tind_ra[1]]);//omgr, orbital frequency in r direction
		orbitpar[5] = (philist[tind_ra[0]] - philist[tind_ra[1]]) / (tlist[tind_ra[0]] - tlist[tind_ra[1]]);//omgphi, orbital frequency in phi direction
	}
	else if (rpi >= 2) {
		orbitpar[3] = 2 * Pi / (tlist[tind_rp[0]] - tlist[tind_rp[1]]);//omgr, orbital frequency in r direction
		orbitpar[5] = (philist[tind_rp[0]] - philist[tind_rp[1]]) / (tlist[tind_rp[0]] - tlist[tind_rp[1]]);//omgphi,orbital frequency in phi direction
	}
	else {
		orbitpar[3] = abs(Pi / (tlist[tind_rp[0]] - tlist[tind_ra[0]]));
		orbitpar[5] = abs((philist[tind_rp[0]] - philist[tind_ra[0]]) / (tlist[tind_rp[0]] - tlist[tind_ra[0]]));
	}


	if (thmaxi >= 2) {
		orbitpar[4] = 2 * Pi / (tlist[tind_thmax[0]] - tlist[tind_thmax[1]]);//omgth, orbital frequency in th direction
	}
	else if (thmini >= 2) {
		orbitpar[4] = 2 * Pi / (tlist[tind_thmin[0]] - tlist[tind_thmin[1]]);//omgth, orbital frequency in th direction
	}
	else {
		orbitpar[4] = abs(Pi / (tlist[tind_thmin[0]] - tlist[tind_thmax[0]]));
	}*/

}


void get_flux(double spin,double E, double L, double Q, double flux[]) {
	//Ref.  PhysRevD.73.064037 
	//double Q_e0 = Q;
	double a = spin, e = current_e, p = current_p, iota = current_iota;
	double p2 = p*p, p3 = pow(p, 3), p4 = pow(p, 4), p5 = pow(p, 5), p6 = pow(p, 6), p7 = pow(p, 7), p8 = pow(p, 8);
	double a2 = a*a, a3 = pow(a, 3), a4 = pow(a, 4), a5 = pow(a, 5), a6 = pow(a, 6), a7 = pow(a, 7), a8 = pow(a, 8);
	double e2 = e*e, e3 = pow(e, 3), e4 = pow(e, 4), e5 = pow(e, 5), e6 = pow(e, 6), e7 = pow(e, 7), e8 = pow(e, 8);

	double g1 = 1 + 73.0 / 24.0*e2 + 37.0 / 96.0*e4, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, ga10, gb10;
	g2 = 73.0 / 12.0 + 823.0 / 24.0* e2 + 949.0 / 32.0*e4 + 491.0 / 192.0 *	e6;
	g3 = 1247.0 / 336.0 + 9181.0 / 672.0 * e2;
	g4 = 4.0 + 1375.0 / 48.0*e2;
	g5 = 44711.0 / 9072.0 + 172157.0 / 2592.0 * e2;
	g6 = 33.0 / 16.0 + 359.0 / 32.0 * e2;
	g7 = 8191.0 / 672.0 + 44531.0 / 336.0 * e2;
	g8 = 3749.0 / 336.0 - 5143.0 / 168.0*  e2;
	g9 = 1.0 + 7.0 / 8.0*e2;
	g10 = 61.0 / 12.0 + 119.0 / 8.0 * e2 + 183.0 / 32.0 * e4;
	g11 = 1247.0 / 336.0 + 425.0 / 336.0 * e2;
	g12 = 4.0 + 97.0 / 8.0 * e2;
	g13 = 44711.0 / 9072.0 + 302893.0 / 6048.0* e2;
	g14 = 33.0 / 16.0 + 95.0 / 16.0* e2;
	g15 = 8191.0 / 672.0 + 48361.0 / 1344.0 * e2;
	g16 = 417.0 / 56.0 - 37241.0 / 672.0 * e2;

	ga10 = 61.0 / 24.0 + 63.0 / 8.0 * e2 + 95.0 / 64.0 *e4;
	gb10 = 61.0 / 8.0 + 91.0 / 4.0 *e2 + 461.0 / 64.0*e4;

	double Eflux2PN, Lflux2PN, Qflux2PN;
	Eflux2PN = -32.0 / 5.0 *massratio / p5 *(1 - e2)*sqrt(1 - e2) *
		(g1 - a / p / sqrt(p) * g2 * cos(iota) - g3 / p + Pi / p / sqrt(p)*g4
			- g5 / p2 + a*a / p2 *g6 - 527.0 / 96.0 *a*a / p2 * pow(sin(iota), 2));//E flux eq.(44) 
	Lflux2PN = -32.0 / 5.0 *massratio / p4 * sqrt(p) *(1 - e2)*sqrt(1 - e2) *
		(g9*cos(iota) + a / p / sqrt(p) * (ga10 - pow(cos(iota), 2)*gb10) - g11 / p*cos(iota)
			+ Pi / p / sqrt(p) * g12 * cos(iota) - g13 / p2*cos(iota) + a*a / p2 * cos(iota) * (g14 - 45.0 / 8.0 * pow(sin(iota), 2)));//Lz flux eq.(45)
																																	   //flux[2] = 2 * Q / L*flux[1];
	Qflux2PN = -64.0 / 5.0 * massratio / p4 * sqrt(p) * sqrt(Q) *sin(iota)*(1 - e2)*sqrt(1 - e2) *
		(g9 - a / p / sqrt(p) *cos(iota) * gb10 - g11 / p + Pi / p / sqrt(p) *g12 - g13 / p2 + a*a / p2*(g14 - 45.0 / 8.0*pow(sin(iota), 2))
			);//eq. (56)

			  //计算e=0时的flux
	e = 0;
	e2 = e*e, e3 = pow(e, 3), e4 = pow(e, 4), e5 = pow(e, 5), e6 = pow(e, 6), e7 = pow(e, 7), e8 = pow(e, 8);

	g1 = 1 + 73.0 / 24.0*e2 + 37.0 / 96.0*e4;
	g2 = 73.0 / 12.0 + 823.0 / 24.0* e2 + 949.0 / 32.0*e4 + 491.0 / 192.0 *	e6;
	g3 = 1247.0 / 336.0 + 9181.0 / 672.0 * e2;
	g4 = 4.0 + 1375.0 / 48.0*e2;
	g5 = 44711.0 / 9072.0 + 172157.0 / 2592.0 * e2;
	g6 = 33.0 / 16.0 + 359.0 / 32.0 * e2;
	g7 = 8191.0 / 672.0 + 44531.0 / 336.0 * e2;
	g8 = 3749.0 / 336.0 - 5143.0 / 168.0*  e2;
	g9 = 1.0 + 7.0 / 8.0*e2;
	g10 = 61.0 / 12.0 + 119.0 / 8.0 * e2 + 183.0 / 32.0 * e4;
	g11 = 1247.0 / 336.0 + 425.0 / 336.0 * e2;
	g12 = 4.0 + 97.0 / 8.0 * e2;
	g13 = 44711.0 / 9072.0 + 302893.0 / 6048.0* e2;
	g14 = 33.0 / 16.0 + 95.0 / 16.0* e2;
	g15 = 8191.0 / 672.0 + 48361.0 / 1344.0 * e2;
	g16 = 417.0 / 56.0 - 37241.0 / 672.0 * e2;

	ga10 = 61.0 / 24.0 + 63.0 / 8.0 * e2 + 95.0 / 64.0 *e4;
	gb10 = 61.0 / 8.0 + 91.0 / 4.0 *e2 + 461.0 / 64.0*e4;

	double Eflux2PN_0, Lflux2PN_0, Qflux2PN_oversqQ0_0, Qfluxsago_0;
	Eflux2PN_0 = -32.0 / 5.0 *massratio / p5 *(1 - e2)*sqrt(1 - e2) *
		(g1 - a / p / sqrt(p) * g2 * cos(iota) - g3 / p + Pi / p / sqrt(p)*g4
			- g5 / p2 + a*a / p2 *g6 - 527.0 / 96.0 *a*a / p2 * pow(sin(iota), 2));//E flux eq.(44) 
	Lflux2PN_0 = -32.0 / 5.0 *massratio / p4 * sqrt(p) *(1 - e2)*sqrt(1 - e2) *
		(g9*cos(iota) + a / p / sqrt(p) * (ga10 - pow(cos(iota), 2)*gb10) - g11 / p*cos(iota)
			+ Pi / p / sqrt(p) * g12 * cos(iota) - g13 / p2*cos(iota) + a*a / p2 * cos(iota) * (g14 - 45.0 / 8.0 * pow(sin(iota), 2)));//Lz flux eq.(45)
																																	   //flux[2] = 2 * Q / L*flux[1];
	Qflux2PN_oversqQ0_0 = -64.0 / 5.0 * massratio / p4 * sqrt(p)  *sin(iota)*(1 - e2)*sqrt(1 - e2) *
		(g9 - a / p / sqrt(p) *cos(iota) * gb10 - g11 / p + Pi / p / sqrt(p) *g12 - g13 / p2 + a*a / p2*(g14 - 45.0 / 8.0*pow(sin(iota), 2))
			);//eq. (56)
	Qfluxsago_0 = -64.0 / 5.0 *massratio / p4 * sqrt(p) *sqrt(Q) *sin(iota) *(1 - e2)*sqrt(1 - e2) *
		(g9 - a / p / sqrt(p) *cos(iota) * gb10 - g11 / p + Pi / p / sqrt(p) *g12 - g13 / p2 + a*a / p2*(211.0 / 96.0 + 847.0 / 96.0*e2)
			);


	double da1 = -10.7420, db1 = 28.5942,
		dc1 = -9.07738, da2 = -1.42836, db2 = 10.7003, dc2 = -33.7090, ca1 = -28.1517, cb1 = 60.9607, cc1 = 40.9998, ca2 = -0.348161,
		cb2 = 2.37258, cc2 = -66.6584, ca3 = -0.715392, cb3 = 3.21593, cc3 = 5.28888,
		ca4 = -7.61034, cb4 = 128.878, cc4 = -475.465, ca5 = 12.2908, cb5 = -113.125, cc5 = 306.119, ca6 = 40.9259, cb6 = -347.271, cc6 = 886.503, ca7 = -25.4831,
		cb7 = 224.227, cc7 = -490.982, ca8 = -9.00634, cb8 = 91.1767, cc8 = -297.002, ca9 = -0.645000, cb9 = -5.13592, cc9 = 47.1982, fa1 = -283.955, fb1 = 736.209,
		fa2 = 483.266, fb2 = -1325.19, fa3 = -219.224, fb3 = 634.499, fa4 = -25.8203, fb4 = 82.0780, fa5 = 301.478, fb5 = -904.161, fa6 = -271.966, fb6 = 827.319;
	double Lfit;
	Lfit = -32.0 / 5.0*massratio / p3 / sqrt(p) *
		(cos(iota) + a / p / sqrt(p)* (61.0 / 24.0 - 61.0 / 8.0* pow(cos(iota), 2)) - 1247.0 / 336.0 / p *cos(iota) + 4 * Pi / p / sqrt(p) *cos(iota) - 44711.0 / 9072.0 / p2 *cos(iota)
			+ a2 / p2  *cos(iota)*  (33.0 / 16.0 - 45.0 / 8.0 *pow(sin(iota), 2))
			+ 1.0 / sqrt(p) / p2*
			(a*(da1 + db1 / sqrt(p) + dc1 / p) + a3 * (da2 + db2 / sqrt(p) + dc2 / p)
				+ cos(iota)*(ca1 + cb1 / sqrt(p) + cc1 / p) + a2 *cos(iota)*(ca2 + cb2 / sqrt(p) + cc2 / p) + a4 *cos(iota)*(ca3 + cb3 / sqrt(p) + cc3 / p)
				+ a *pow(cos(iota), 2)*(ca4 + cb4 / sqrt(p) + cc4 / p) + a3 *pow(cos(iota), 2)*(ca5 + cb5 / sqrt(p) + cc5 / p) + a2 *pow(cos(iota), 3)*(ca6 + cb6 / sqrt(p) + cc6 / p)
				+ a4 *pow(cos(iota), 3)*(ca7 + cb7 / sqrt(p) + cc7 / p) + a3 *pow(cos(iota), 4)*(ca8 + cb8 / sqrt(p) + cc8 / p) + a4 *pow(cos(iota), 5)*(ca9 + cb9 / sqrt(p) + cc9 / p))
			+ 1.0 / sqrt(p) / p3 * a *cos(iota)*
			(fa1 + fb1 / sqrt(p) + a*(fa2 + fb2 / sqrt(p)) + a2*(fa3 + fb3 / sqrt(p)) + pow(cos(iota), 2)*(fa4 + fb4 / sqrt(p))
				+ a *pow(cos(iota), 2)*(fa5 + fb5 / sqrt(p)) + a2 *pow(cos(iota), 2)*(fa6 + fb6 / sqrt(p)))
			);//eq.(57)

	double c10a = -0.0309341, c10b = -22.2416, c10c = 7.55265, c11a = -3.33476, c11b = 22.7013, c11c = -12.4700,
		fa7 = -162.268, fb7 = 247.168, fa8 = 152.125, fb8 = -182.165, fa9 = 184.465, fb9 = -267.553, fa10 = -188.132, fb10 = 254.067;
	double iotafit_sqQe0 = 32.0 / 5.0 * massratio * a *pow(sin(iota), 2) / p5 * (
		61.0 / 24.0 + 1.0 / p * (da1 + db1 / sqrt(p) + dc1 / p)
		+ a2 / p*(da2 + db2 / sqrt(p) + dc2 / p) + a*cos(iota) / sqrt(p) * (c10a + c10b / p + c10c / p / sqrt(p))
		+ a2*pow(cos(iota), 2) / p * (c11a + c11b / sqrt(p) + c11c / p)
		+ 1.0 / p2 / sqrt(p) *a3 * cos(iota) * (fa7 + fb7 / sqrt(p)
			+ a* (fa8 + fb8 / sqrt(p)) + pow(cos(iota), 2) *(fa9 + fb9 / sqrt(p)) + a * pow(cos(iota), 2) * (fa10 + fb10 / sqrt(p)))
		);
	//printf("iotadot=%f\n", iotafit_sqQe0/sqrt(Q_e0));
	double Qfit = 2 * tan(iota)*sqrt(Q)*(Lfit + iotafit_sqQe0 / pow(sin(iota), 2));
	e = current_e;
	e2 = e*e, e3 = pow(e, 3), e4 = pow(e, 4), e5 = pow(e, 5), e6 = pow(e, 6), e7 = pow(e, 7), e8 = pow(e, 8);
	double Lfluxmod = Lflux2PN + pow(1 - e*e, 1.5) * (-Lflux2PN_0 + Lfit);//eq.(59)
	double Qfluxmod = Qflux2PN + pow(1 - e*e, 1.5) * (-sqrt(Q) * Qflux2PN_oversqQ0_0 + Qfit);
	double N1, N4, N5;
	N1 = E*p4 + a*a*E*p2 - 2 * a*(L - a*E)*p;//eq. (10)
	N4 = (2.0 * p - p2)*L - 2.0 * a*E*p;//eq.(15)
	N5 = (2 * p - p2 - a2) / 2.0;//eq.(15)
	double Efluxmod = Eflux2PN - pow(1 - e*e, 1.5) * (Eflux2PN_0 + N4 / N1*Lfit + N5 / N1*Qfit);//eq.(20)
	flux[0] = Efluxmod;
	flux[1] = Lfluxmod;
	flux[2] = Qfluxmod;
}

void dcardboy(double spin, double boy[], double car[], double dcar_dboy[][4]) {//从Boyer-Lindquist坐标换到Cartesian坐标
	car[0] = boy[0];
	car[1] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) *cos(boy[3]);
	car[2] = sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) *sin(boy[3]);
	car[3] = boy[1] * cos(boy[2]);
	dcar_dboy[0][0] = 1;
	dcar_dboy[0][1] = 0;
	dcar_dboy[0][2] = 0;
	dcar_dboy[0][3] = 0;
	dcar_dboy[1][0] = 0;
	dcar_dboy[2][0] = 0;
	dcar_dboy[3][0] = 0;
	dcar_dboy[1][1] = boy[1] / sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * cos(boy[3]);//dx/dr
																							 //好像不太方便全用boy[]
	double r = boy[1], th = boy[2], phi = boy[3];
	double a = spin;
	dcar_dboy[1][2] = sqrt(r*r + spin*spin)*cos(th)*cos(phi);//dx/dth
	dcar_dboy[1][3] = -sqrt(r*r + spin*spin)*sin(th)*sin(phi);// dx/dphi
	dcar_dboy[2][1] = boy[1] / sqrt(boy[1] * boy[1] + spin*spin) * sin(boy[2]) * sin(boy[3]);//dy/dr
	dcar_dboy[2][2] = sqrt(r*r + a*a)*cos(th)*sin(phi);// dy/dth
	dcar_dboy[2][3] = sqrt(r*r + a*a)*sin(th)*cos(phi);// dy/dphi
	dcar_dboy[3][1] = cos(th);//dz/dr
	dcar_dboy[3][2] = -r*sin(th);//dz/dth
	dcar_dboy[3][3] = 0;
}
void dboydcar(double spin, double car[], double boy[], double dboy_dcar[][4]) {//从Cartesian坐标换到Boyer-Lindquist坐标
	double x = car[1],y = car[2],z = car[3];
	double a = spin;
	double r = sqrt(sqrt((-a *a + x *x + y *y + z *z)*(-a *a + x *x + y *y + z *z) + 4 * a *a * z*z) - a*a + x *x + y*y + z*z) / sqrt(2);
	double th;
	if (a == 0) {
		th = acos(z / sqrt(x*x + y*y + z*z));
	}
	else {
		double sqcth;
		if (abs(z*z / (x*x + y*y + a*a)) > 1e-10) {
			sqcth = (1 - x*x / a / a - y*y / a / a - z*z / a / a + sqrt(4 * a*a*z*z + (x*x + y*y + z*z - a*a)*(x*x + y*y + z*z - a*a)) / a / a) / 2;
			if (abs(sqcth) < eps) th = Piby2-sqrt(sqcth);
			else {
				if (z > 0) {
					double cth = sqrt(abs(sqcth));
					th = acos(cth);
				}
				else {
					double cth = -sqrt(abs(sqcth));
					th = acos(cth);
				}
			}
		}
		else {// z*z/(x*x+y*y+a*a) 趋于0 的近似
			double cth = z/ a *sqrt((a*a + z*z / 4) / (x*x + y*y - a*a));
			sqcth = cth*cth;
			th = acos(cth);
		}
	}
	double phi;
	if (x == 0) {
		if (y>0) phi = Piby2;
		else phi = -Piby2;
	}
	else {
		phi = atan(y / x);
		if (x < 0) phi += Pi;//atan的值域在(-pi/2, pi/2)
	}

	boy[0] = car[0];boy[1] = r;boy[2] = th;boy[3] = phi;
	dboy_dcar[0][0] = 1;
	dboy_dcar[0][1] = 0;
	dboy_dcar[0][2] = 0;
	dboy_dcar[0][3] = 0;
	dboy_dcar[1][0] = 0;
	dboy_dcar[2][0] = 0;
	dboy_dcar[3][0] = 0;
	dboy_dcar[1][1] = r*r*r*x/( r*r*r*r+a*a*z*z );//dr/dx
	dboy_dcar[1][2] = r*r*r*y / (r*r*r*r + a*a*z*z);//dr/dy
	dboy_dcar[1][3] = r*(a*a + r*r)*z / (r*r*r*r + a*a*z*z);//dr/dz
	
	if (abs(th - Piby2) < eps) { //th趋于pi/2的极限
		dboy_dcar[2][1] = 0;
		dboy_dcar[2][2] = 0;
		dboy_dcar[2][3] = - 1 / r;
	}
	else if (abs(th - 0) < eps) { //th趋于0的极限
		dboy_dcar[2][1] = cos(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][2] = sin(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][3] = 0;
	}
	else {
		dboy_dcar[2][1] = x*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][2] = y*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][3] = -z *sin(th)*sin(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));

	}

	
	dboy_dcar[3][1] = -y / (x*x + y*y);
	dboy_dcar[3][2] = x / (x*x + y*y);
	//dphi/dx
	//dphi/dy
	dboy_dcar[3][3] = 0;//dphi/dz
	
}

void dboydcar_boyknown(double spin, double car[], double boy[], double dboy_dcar[][4]) {//已知两套坐标，快速一点算导数
	double x = car[1], y = car[2], z = car[3];
	double a = spin;
	double r = boy[1], th = boy[2], phi = boy[3];

	dboy_dcar[0][0] = 1;
	dboy_dcar[0][1] = 0;
	dboy_dcar[0][2] = 0;
	dboy_dcar[0][3] = 0;
	dboy_dcar[1][0] = 0;
	dboy_dcar[2][0] = 0;
	dboy_dcar[3][0] = 0;
	dboy_dcar[1][1] = r*r*r*x / (r*r*r*r + a*a*z*z);//dr/dx
	dboy_dcar[1][2] = r*r*r*y / (r*r*r*r + a*a*z*z);//dr/dy
	dboy_dcar[1][3] = r*(a*a + r*r)*z / (r*r*r*r + a*a*z*z);//dr/dz

	if (abs(th - Piby2) < eps) { //th趋于pi/2的极限
		dboy_dcar[2][1] = 0;
		dboy_dcar[2][2] = 0;
		dboy_dcar[2][3] = -1 / r;
	}
	else if (abs(th - 0) < eps) { //th趋于0的极限
		dboy_dcar[2][1] = cos(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][2] = sin(phi) / sqrt(r*r + a*a);
		dboy_dcar[2][3] = 0;
	}
	else {
		dboy_dcar[2][1] = x*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][2] = y*cos(th)*cos(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));
		dboy_dcar[2][3] = -z *sin(th)*sin(th) / (a*a*cos(th)*cos(th)*cos(th)*sin(th) + z*z*tan(th));

	}

	if (abs(x - 0) < eps) {//x趋于0的极限
		dboy_dcar[3][1] = -1 / y;
		dboy_dcar[3][2] = 0;
	}
	else if (abs(y - 0)<eps) {//y趋于0的极限
		dboy_dcar[3][1] = 0;
		dboy_dcar[3][2] = 1 / x;
	}
	else {
		dboy_dcar[3][1] = -y / (x*x + y*y);
		dboy_dcar[3][2] = x / (x*x + y*y);
	}
	//dphi/dx
	//dphi/dy
	dboy_dcar[3][3] = 0;//dphi/dz

}

void highderiv(double spin, double defpar[], double boy[], double boy_taudot[], double x[][4]) {//坐标对t的s阶导数存在x[s]里
	//写的时候注意保留x[][]为Cartesian坐标下对t的0阶至7阶导数，Gamma为Cartesian坐标下的Gamma

	double dcar_dboy[4][4] = { 0 }, dboy_dcar[4][4] = { 0 };
	double Gamma_boy[4][4][4],Gamma[4][4][4];
	double r = boy[1], th = boy[2];
	dcardboy(spin, boy, x[0], dcar_dboy);//算出cartesian坐标下的0阶导数
	dboydcar_boyknown(spin, x[0], boy, dboy_dcar);
	Christoffel(spin, defpar, r, th, Gamma_boy);
	int al, bt, gm, mu, sg, ro;//指标alpha, beta, gamma, mu, sigma, rho
	
	//转换Christoffel到Cartesian 坐标 {\Gamma^{(Cart)} }^\alpha_{\beta \gamma} = \Gamma^{(Boyer)} ^\mu_{\sigma\rho} \frac{\partial x^{(cart)}^\alpha }{x^{(boyer)} ^\mu } \frac{\partial x^{(boyer)}^\sigma }{x^{(cart)} ^\beta } \frac{\partial x^{(boyer)}^\rho }{x^{(cart)} ^\gamma }
	for (al = 0;al < 4;al++) {
		for (bt = 0;bt < 4;bt++) {
			for (gm = 0;gm < 4;gm++) {
				Gamma[al][bt][gm] = 0;
				for (mu = 0;mu < 4;mu++) {
					for (sg = 0;sg < 4;sg++) {
						for (ro = 0;ro < 4;ro++) {
							Gamma[al][bt][gm] += Gamma_boy[mu][sg][ro] * dcar_dboy[al][mu] * dboy_dcar[sg][bt] * dboy_dcar[ro][gm];
						}
					}
				}
			}
		}
	}

	double car_taudot[4];//速度也换成Cartesian坐标
	for (al = 0;al < 4;al++) {
		car_taudot[al] = 0;
		for (mu = 0;mu < 4;mu++) {
			car_taudot[al] += boy_taudot[mu] * dcar_dboy[al][mu];
		}
	}

	//对t的1阶导数
	x[1][0] = 1;x[1][1] = car_taudot[1] / car_taudot[0]; x[1][2] = car_taudot[2] / car_taudot[0]; x[1][3] = car_taudot[3] / car_taudot[0];

	//对t的2阶导数
	x[2][0] = 0;
	int i,j,k,l;//遍历1到3的指标i
	int nu;
	for (i = 1;i < 4;i++) {
		x[2][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[2][i] += (-Gamma[i][mu][nu] * x[1][mu] * x[1][nu] + Gamma[0][mu][nu] * x[1][mu] * x[1][nu] * x[1][i]);
			}
		}
	}

	//后面需要算Gamma的导数了
	/*
	double Gamma_dt[5][4][4][4] = { 0 };//Gamma_dt存的是对t的各阶导数
	double Gamma_der1[4][4][4][4] = { 0 };//Gamma_ders存的是Gamma对坐标的s阶导数
	double derstep = 5e-4;//求导的eps
	*/
	//Gamma对t的一阶导数
	/*
	double xtemp[4];//临时存放坐标
	for (i = 1;i < 4;i++) {
		for (int ind = 0;ind < 4;ind++) xtemp[ind] = x[0][ind];
		xtemp[i] += derstep / 2;
		double Gamma_ip[4][4][4] = { 0 };
		Gamma_car(spin, defpar, xtemp, Gamma_ip);
		xtemp[i] -= derstep ;
		double Gamma_im[4][4][4] = { 0 };
		Gamma_car(spin, defpar, xtemp, Gamma_im);
		for (al = 0;al < 4;al += 1) {
			for (bt = 0;bt < 4;bt += 1) {
				for (gm = 0;gm < 4;gm += 1) {
					Gamma_der1[i][al][bt][gm] = (Gamma_ip[al][bt][gm] - Gamma_im[al][bt][gm]) / derstep;//计算对坐标的导数
					Gamma_dt[1][al][bt][gm] += Gamma_der1[i][al][bt][gm] * x[1][i];//对i求和 得到对t的导数
				}
			}
		}

	}
	*/
	//Gamma对t的二阶导数
	/*
	double Gamma_der2[4][4][4][4][4] = { 0 };
	for (i = 1;i < 4;i++) {
		for (j = 1;j < 4;j++) {
			for (int ind = 0;ind < 4;ind++) xtemp[ind] = x[0][ind];
			double Gamma_ipjp[4][4][4], Gamma_ipjm[4][4][4], Gamma_imjp[4][4][4], Gamma_imjm[4][4][4];
			xtemp[i] += derstep / 2;xtemp[j] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];
			xtemp[i] += derstep / 2;xtemp[j] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];
			xtemp[i] -= derstep / 2;xtemp[j] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];
			xtemp[i] -= derstep / 2;xtemp[j] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];
			for (al = 0;al < 4;al += 1) {
				for (bt = 0;bt < 4;bt += 1) {
					for (gm = 0;gm < 4;gm += 1) {
						Gamma_der2[i][j][al][bt][gm] = (Gamma_ipjp[al][bt][gm] + Gamma_imjm[al][bt][gm] - Gamma_imjp[al][bt][gm] - Gamma_ipjm[al][bt][gm]) / pow(derstep, 2);//计算对坐标的导数
						Gamma_dt[2][al][bt][gm] += Gamma_der2[i][j][al][bt][gm] * x[1][i] * x[1][j];//对ij求和 得到对t的导数
						if (j == 1) {//不对j求和
							Gamma_dt[2][al][bt][gm] += Gamma_der1[i][al][bt][gm] * x[2][i];
						}
					}
				}
			}
		}
	}
	*/
	//Gamma对t的三阶导数
	/*
	double Gamma_der3[4][4][4][4][4][4] = { 0 };
	for (i = 1;i < 4;i++) {
		for (j = 1;j < 4;j++) {
			for (k = 1;k < 4;k++) {
				for (int ind = 0;ind < 4;ind++) xtemp[ind] = x[0][ind];
				double Gamma_ipjpkp[4][4][4] = { 0 }, Gamma_ipjpkm[4][4][4] = { 0 }, Gamma_ipjmkp[4][4][4] = { 0 }, Gamma_ipjmkm[4][4][4] = { 0 };
				double Gamma_imjpkp[4][4][4] = { 0 }, Gamma_imjpkm[4][4][4] = { 0 }, Gamma_imjmkp[4][4][4] = { 0 }, Gamma_imjmkm[4][4][4] = { 0 };
				xtemp[i] += derstep / 2;xtemp[j] += derstep / 2;xtemp[k] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjpkp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] += derstep / 2;xtemp[j] += derstep / 2;xtemp[k] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjpkm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] += derstep / 2;xtemp[j] -= derstep / 2;xtemp[k] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjmkp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] += derstep / 2;xtemp[j] -= derstep / 2;xtemp[k] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_ipjmkm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] -= derstep / 2;xtemp[j] += derstep / 2;xtemp[k] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjpkp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] -= derstep / 2;xtemp[j] += derstep / 2;xtemp[k] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjpkm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] -= derstep / 2;xtemp[j] -= derstep / 2;xtemp[k] += derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjmkp);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				xtemp[i] -= derstep / 2;xtemp[j] -= derstep / 2;xtemp[k] -= derstep / 2;Gamma_car(spin, defpar, xtemp, Gamma_imjmkm);xtemp[i] = x[0][i];xtemp[j] = x[0][j];xtemp[k] = x[0][k];
				for (al = 0;al < 4;al += 1) {
					for (bt = 0;bt < 4;bt += 1) {
						for (gm = 0;gm < 4;gm += 1) {
							Gamma_der3[i][j][k][al][bt][gm] = (Gamma_ipjpkp[al][bt][gm] - Gamma_ipjpkm[al][bt][gm] - Gamma_ipjmkp[al][bt][gm] + Gamma_ipjmkm[al][bt][gm]
								-Gamma_imjpkp[al][bt][gm] + Gamma_imjpkm[al][bt][gm] + Gamma_imjmkp[al][bt][gm] - Gamma_imjmkm[al][bt][gm]) / pow(derstep, 3);//计算对坐标的导数
							Gamma_dt[3][al][bt][gm] += Gamma_der3[i][j][k][al][bt][gm] * x[1][i] * x[1][j] * x[1][k];//对ijk求和 得到对t的导数
							if (j == 1&&k==1) {//不对ij求和
								Gamma_dt[3][al][bt][gm] += Gamma_der1[i][al][bt][gm] * x[3][i];
							}
							if (k == 1) {//不对k求和
								Gamma_dt[3][al][bt][gm] += Gamma_der2[i][j][al][bt][gm] * x[1][i] * x[2][j];
							}
						}
					}
				}
			}
		}
	}
	*/

	//对t的3阶导数
	x[3][0] = 0;
	for (i = 1;i < 4;i++) {
		x[3][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[3][i] += (-Gamma[i][mu][nu] * (x[2][mu] * x[1][nu] + x[1][mu] * x[2][nu])
					//- Gamma_dt[1][i][mu][nu] * x[1][mu] * x[1][nu]
					+ Gamma[0][mu][nu] * (x[2][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[2][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[2][i]))
					/*+ Gamma_dt[1][0][mu][nu]* x[1][mu] * x[1][nu] * x[1][i]*/ ;
			}
		}
	}

	//对t的4阶导数
	x[4][0] = 0;
	for (i = 1;i < 4;i++) {
		x[4][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[4][i] += (-Gamma[i][mu][nu] * (x[3][mu] * x[1][nu] + 2 * x[2][mu] * x[2][nu] + x[1][mu] * x[3][nu])
					+ Gamma[0][mu][nu] * (x[3][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[3][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[3][i]
						+ 2 * x[2][mu] * x[2][nu] * x[1][i] + 2 * x[1][mu] * x[2][nu] * x[2][i] + 2 * x[2][mu] * x[1][nu] * x[2][i]));
			}
		}
	}

	//对t的5阶导数
	x[5][0] = 0;
	for (i = 1;i < 4;i++) {
		x[5][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[5][i] += (-Gamma[i][mu][nu] * (x[4][mu] * x[1][nu] + 3 * x[3][mu] * x[2][nu] + 3* x[2][mu]*x[3][nu] + x[1][mu] * x[4][nu])
					+ Gamma[0][mu][nu] * (x[4][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[4][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[4][i]
						+ 3 * x[3][mu] * x[2][nu] * x[1][i] + 3 * x[2][mu] * x[3][nu] * x[1][i] + 3 * x[1][mu] * x[3][nu] * x[2][i] + 3 * x[1][mu] * x[2][nu] * x[3][i] 
						+ 3* x[2][mu]*x[1][nu]*x[3][i] +3* x[3][mu]*x[1][nu]*x[2][i] + 6 * x[2][mu]*x[2][nu]*x[2][i] ));
			}
		}
	}

	//对t的6阶导数
	x[6][0] = 0;
	for (i = 1;i < 4;i++) {
		x[6][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[6][i] += (-Gamma[i][mu][nu] * (x[5][mu] * x[1][nu] + 4 * x[4][mu] * x[2][nu] + 6 * x[3][mu]*x[3][nu] + 4 * x[2][mu] * x[4][nu] + x[1][mu] * x[5][nu])
					+ Gamma[0][mu][nu] * (x[5][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[5][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[5][i]
						+ 4 * x[4][mu] * x[2][nu] * x[1][i] + 4 * x[2][mu] * x[4][nu] * x[1][i] + 4 * x[1][mu] * x[4][nu] * x[2][i]
						+ 4 * x[1][mu] * x[2][nu] * x[4][i] + 4 * x[2][mu] * x[1][nu] * x[4][i] + 4 * x[4][mu] * x[1][nu] * x[2][i]
						+ 6 * x[3][mu] * x[3][nu] * x[1][i] + 6 * x[3][mu] * x[1][nu] * x[3][i] + 6 * x[1][mu] * x[3][nu] * x[3][i]
						+ 12* x[2][mu] * x[2][nu] * x[3][i] + 12* x[2][mu] * x[3][nu] * x[2][i] + 12* x[3][mu] * x[2][nu] * x[2][i] ));
			}
		}
	}

	//对t的7阶导数
	x[7][0] = 0;
	for (i = 1;i < 4;i++) {
		x[7][i] = 0;
		for (mu = 0;mu < 4;mu++) {
			for (nu = 0;nu < 4;nu++) {
				x[7][i] += (-Gamma[i][mu][nu] * (x[5][mu] * x[1][nu] + 4 * x[4][mu] * x[2][nu] + 6 * x[3][mu] * x[3][nu] + 4 * x[2][mu] * x[4][nu] + x[1][mu] * x[5][nu])
					+ Gamma[0][mu][nu] * (x[5][mu] * x[1][nu] * x[1][i] + x[1][mu] * x[5][nu] * x[1][i] + x[1][mu] * x[1][nu] * x[5][i]
						+ 4 * x[4][mu] * x[2][nu] * x[1][i] + 4 * x[2][mu] * x[4][nu] * x[1][i] + 4 * x[1][mu] * x[4][nu] * x[2][i]
						+ 4 * x[1][mu] * x[2][nu] * x[4][i] + 4 * x[2][mu] * x[1][nu] * x[4][i] + 4 * x[4][mu] * x[1][nu] * x[2][i]
						+ 6 * x[3][mu] * x[3][nu] * x[1][i] + 6 * x[3][mu] * x[1][nu] * x[3][i] + 6 * x[1][mu] * x[3][nu] * x[3][i]
						+ 12 * x[2][mu] * x[2][nu] * x[3][i] + 12 * x[2][mu] * x[3][nu] * x[2][i] + 12 * x[3][mu] * x[2][nu] * x[2][i]));
			}
		}
	}
}

void Gamma_highderiv(double spin,double defpar[],double x[][4],double Gamma_dt[][4][4][4]) {//Gamma_dt存的是对t的各阶导数
	double Gamma_der1[4][4][4][4];//Gamma_ders存的是Gamma对坐标的s阶导数


}

void Gamma_car(double spin, double defpar[],double car[],double Gamma[][4][4]) {//Cartesian坐标下的Gamma
	double dcar_dboy[4][4] = { 0 }, dboy_dcar[4][4] = { 0 };
	double Gamma_boy[4][4][4] = { 0 };
	double boy[4] = { 0 }, cartemp[4] = { 0 };
	dboydcar(spin, car, boy, dboy_dcar);
	dcardboy(spin, boy, cartemp, dcar_dboy);
	double r = boy[1], th = boy[2];
	Christoffel(spin, defpar, r, th, Gamma_boy);
	int al, bt, gm, mu, sg, ro;//指标alpha, beta, gamma, mu, sigma, rho

							   //转换Christoffel到Cartesian 坐标 {\Gamma^{(Cart)} }^\alpha_{\beta \gamma} = \Gamma^{(Boyer)} ^\mu_{\sigma\rho} \frac{\partial x^{(cart)}^\alpha }{x^{(boyer)} ^\mu } \frac{\partial x^{(boyer)}^\sigma }{x^{(cart)} ^\beta } \frac{\partial x^{(boyer)}^\rho }{x^{(cart)} ^\gamma }
	for (al = 0;al < 4;al++) {
		for (bt = 0;bt < 4;bt++) {
			for (gm = 0;gm < 4;gm++) {
				Gamma[al][bt][gm] = 0;
				for (mu = 0;mu < 4;mu++) {
					for (sg = 0;sg < 4;sg++) {
						for (ro = 0;ro < 4;ro++) {
							Gamma[al][bt][gm] += Gamma_boy[mu][sg][ro] * dcar_dboy[al][mu] * dboy_dcar[sg][bt] * dboy_dcar[ro][gm];
						}
					}
				}
			}
		}
	}
}

void radacc(double spin, double x[][4], double acc[]) {//注意，乘上质量比才是真·加速度
	double lev[4][4][4] = { 0 };
	lev[1][2][3] = 1; lev[2][3][1] = 1; lev[3][1][2] = 1;
	lev[3][2][1] = -1; lev[2][1][3] = -1; lev[1][3][2] = -1;//levi-civita记号

	int i, j, k, p, q;//求和的指标
	double I_5dot[4][4] = { 0 }, J_5dot[4][4] = { 0 }, J_6dot[4][4] = { 0 };//0分量放着占位

	//I_5dot[1][1];
	//I_5dot[2][2];
	//I_5dot[3][3];
	//I_5dot[1][2];
	//I_5dot[2][3];
	//I_5dot[1][3];
	for (i = 1;i < 4;i++) {
		for (j = i;j < 4;j++) {
			I_5dot[i][j] = x[5][i] * x[0][j] + 5 * x[4][i] * x[1][j] + 10 * x[3][i] * x[2][j] + 5 * x[2][i] * x[3][j] + x[0][i] * x[5][j];
		}
	}
	I_5dot[2][1] = I_5dot[1][2];
	I_5dot[3][2] = I_5dot[2][3];
	I_5dot[3][1] = I_5dot[1][3];

	//算J的5,6次导数
	int k1, m1, k2, m2;
	for (i = 1;i < 4;i++) {
		for (j = 1;j < 4;j++) {

			if (j == 1) {
				k1 = 2;m1 = 3;
				k2 = 3;m2 = 2;
			}
			else if (j == 2) {
				k1 = 3;m1 = 1;
				k2 = 1;m2 = 3;
			}
			else {
				J_5dot[i][j] = -3.0 / 2.0*spin*x[5][i];
				J_6dot[i][j] = -3.0 / 2.0*spin*x[6][i];//代替delta_j3
				k1 = 1;m1 = 2;
				k2 = 2;m2 = 1;
			}//用来代替levi-civita

			J_5dot[i][j] += (x[5][i] * x[0][k1] * x[1][m1] + x[0][i] * x[5][k1] * x[1][m1] + x[0][i] * x[0][k1] * x[6][m1]
				+ 5 * x[4][i] * x[1][k1] * x[1][m1] + 5 * x[1][i] * x[4][k1] * x[1][m1] + 5 * x[0][i] * x[4][k1] * x[2][m1] + 5 * x[4][i] * x[0][k1] * x[2][m1] + 5 * x[1][i] * x[0][k1] * x[5][m1] + 5 * x[0][i] * x[1][k1] * x[5][m1]
				+ 10 * x[3][i] * x[2][k1] * x[1][m1] + 10 * x[2][i] * x[3][k1] * x[1][m1] + 10 * x[3][i] * x[0][k1] * x[3][m1] + 10 * x[0][i] * x[3][k1] * x[3][m1] + 10 * x[2][i] * x[0][k1] * x[4][m1] + 10 * x[0][i] * x[2][k1] * x[4][m1]
				+ 20 * x[3][i] * x[1][k1] * x[2][m1] + 20 * x[1][i] * x[3][k1] * x[2][m1] + 20 * x[1][i] * x[1][k1] * x[4][m1]
				+ 30 * x[2][i] * x[2][k1] * x[2][m1] + 30 * x[2][i] * x[1][k1] * x[3][m1] + 30 * x[1][i] * x[2][k1] * x[3][m1])
				- (x[5][i] * x[0][k2] * x[1][m2] + x[0][i] * x[5][k2] * x[1][m2] + x[0][i] * x[0][k2] * x[6][m2]
					+ 5 * x[4][i] * x[1][k2] * x[1][m2] + 5 * x[1][i] * x[4][k2] * x[1][m2] + 5 * x[0][i] * x[4][k2] * x[2][m2] + 5 * x[4][i] * x[0][k2] * x[2][m2] + 5 * x[1][i] * x[0][k2] * x[5][m2] + 5 * x[0][i] * x[1][k2] * x[5][m2]
					+ 10 * x[3][i] * x[2][k2] * x[1][m2] + 10 * x[2][i] * x[3][k2] * x[1][m2] + 10 * x[3][i] * x[0][k2] * x[3][m2] + 10 * x[0][i] * x[3][k2] * x[3][m2] + 10 * x[2][i] * x[0][k2] * x[4][m2] + 10 * x[0][i] * x[2][k2] * x[4][m2]
					+ 20 * x[3][i] * x[1][k2] * x[2][m2] + 20 * x[1][i] * x[3][k2] * x[2][m2] + 20 * x[1][i] * x[1][k2] * x[4][m2]
					+ 30 * x[2][i] * x[2][k2] * x[2][m2] + 30 * x[2][i] * x[1][k2] * x[3][m2] + 30 * x[1][i] * x[2][k2] * x[3][m2]);
			J_6dot[i][j]+= (x[6][i] * x[0][k1] * x[1][m1] + x[0][i] * x[6][k1] * x[1][m1] + x[0][i] * x[0][k1] * x[7][m1]
				+ 6 * x[6][i] * x[1][k1] * x[1][m1] + 6 * x[1][i] * x[6][k1] * x[1][m1] + 6 * x[0][i] * x[5][k1] * x[2][m1] + 6 * x[5][i] * x[0][k1] * x[2][m1] + 6 * x[1][i] * x[0][k1] * x[6][m1] + 6 * x[0][i] * x[1][k1] * x[6][m1]
				+ 15 * x[4][i] * x[2][k1] * x[1][m1] + 15 * x[2][i] * x[4][k1] * x[1][m1] + 15 * x[4][i] * x[0][k1] * x[3][m1] + 15 * x[0][i] * x[4][k1] * x[3][m1] + 15 * x[2][i] * x[0][k1] * x[5][m1] + 15 * x[0][i] * x[2][k1] * x[5][m1]
				+ 20 * x[3][i] * x[3][k1] * x[1][m1] + 20 * x[3][i] * x[0][k1] * x[4][m1] + 20 * x[1][i] * x[3][k1] * x[4][m1]
				+ 30 * x[4][i] * x[1][k1] * x[2][m1] + 30 * x[1][i] * x[4][k1] * x[2][m1] + 30 * x[1][i] * x[1][k1] * x[5][m1]
				+ 60 * x[3][i] * x[2][k1] * x[2][m1] + 60 * x[2][i] * x[3][k1] * x[2][m1] + 60 * x[3][i] * x[1][k1] * x[3][m1]
				+ 60 * x[1][i] * x[3][k1] * x[3][m1] + 60 * x[2][i] * x[1][k1] * x[4][m1] + 60 * x[1][i] * x[2][k1] * x[4][m1]
				+ 90 * x[2][i] * x[2][k1] * x[3][m1])
				- (x[6][i] * x[0][k2] * x[1][m2] + x[0][i] * x[6][k2] * x[1][m2] + x[0][i] * x[0][k2] * x[7][m2]
					+ 6 * x[6][i] * x[1][k2] * x[1][m2] + 6 * x[1][i] * x[6][k2] * x[1][m2] + 6 * x[0][i] * x[5][k2] * x[2][m2] + 6 * x[5][i] * x[0][k2] * x[2][m2] + 6 * x[1][i] * x[0][k2] * x[6][m2] + 6 * x[0][i] * x[1][k2] * x[6][m2]
					+ 15 * x[4][i] * x[2][k2] * x[1][m2] + 15 * x[2][i] * x[4][k2] * x[1][m2] + 15 * x[4][i] * x[0][k2] * x[3][m2] + 15 * x[0][i] * x[4][k2] * x[3][m2] + 15 * x[2][i] * x[0][k2] * x[5][m2] + 15 * x[0][i] * x[2][k2] * x[5][m2]
					+ 20 * x[3][i] * x[3][k2] * x[1][m2] + 20 * x[3][i] * x[0][k2] * x[4][m2] + 20 * x[1][i] * x[3][k2] * x[4][m2]
					+ 30 * x[4][i] * x[1][k2] * x[2][m2] + 30 * x[1][i] * x[4][k2] * x[2][m2] + 30 * x[1][i] * x[1][k2] * x[5][m2]
					+ 60 * x[3][i] * x[2][k2] * x[2][m2] + 60 * x[2][i] * x[3][k2] * x[2][m2] + 60 * x[3][i] * x[1][k2] * x[3][m2]
					+ 60 * x[1][i] * x[3][k2] * x[3][m2] + 60 * x[2][i] * x[1][k2] * x[4][m2] + 60 * x[1][i] * x[2][k2] * x[4][m2]
					+ 90 * x[2][i] * x[2][k2] * x[3][m2]);
		}
	}

	//Symmetric-trace-free
	for (i = 1;i < 4;i++) {
		for (j = 1;j < 4;j++) {
			J_5dot[i][j] = (J_5dot[i][j] + J_5dot[j][i]) / 2.0;
			J_6dot[i][j] = (J_5dot[i][j] + J_6dot[j][i]) / 2.0;
			if (i == j) {
				I_5dot[i][j] += -(I_5dot[1][1] + I_5dot[2][2] + I_5dot[3][3]) / 3.0;
				J_5dot[i][j] += -(J_5dot[1][1] + J_5dot[2][2] + J_5dot[3][3]) / 3.0;
				J_6dot[i][j] += -(J_6dot[1][1] + J_6dot[2][2] + J_6dot[3][3]) / 3.0;
			}
		}
	}


	//算加速度
	acc[0] = 0;
	for (j = 1;j < 4;j++) {
		acc[j] = 8.0 / 15.0*spin*J_5dot[3][j];
		for (k = 1;k < 4;k++) {
			acc[j] += -2.0 / 5.0*I_5dot[j][k] * x[0][k];
			for (p = 1;p < 4;p++) {
				for (q = 1;q < 4;q++) {
					acc[j] += 16.0 / 45.0 * lev[j][p][q] * J_6dot[p][k] * x[0][q] * x[0][k] + 32.0 / 45.0 * lev[j][p][q] * J_5dot[p][k] * x[0][k] * x[1][q]
						+16.0/45.0 * lev[p][q][j] * J_5dot[k][p] * x[0][q] *x[1][k] - 16.0/45.0 * lev[p][q][k] * J_5dot[j][p] * x[0][q] * x[1][k];
				}
			}
		}
	}
}
//一个以前用的RK4算法
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
//		metric(spin, defpar, r, th, g);
//		Christoffel(spin, defpar, r, th, Gamma);
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
//		metric(spin, defpar, varnew1[1], varnew1[2], g);
//		Christoffel(spin, defpar, varnew1[1], varnew1[2], Gamma);
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
//
//		/*********************CALCULATE k3********************/
//		metric(spin, defpar, varnew2[1], varnew2[2], g);
//		Christoffel(spin, defpar, varnew2[1], varnew2[2], Gamma);
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
//		metric(spin, defpar, varnew3[1], varnew3[2], g);
//		Christoffel(spin, defpar, varnew3[1], varnew3[2], Gamma);
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
