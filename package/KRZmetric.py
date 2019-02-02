
# coding: utf-8
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import numpy as np
import scipy
from numpy import sin,cos,sqrt
#z1就是r，z2就是theta

def metric_KRZ(spin,   d1,   z1,   z2):
    r = z1;
    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    th = z2;

    sinth = sin(th);
    costh = cos(th);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    fourthspin = sqspin*sqspin;
    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;
    k23 = sqspin / sqr0;

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber)+(a20*cuber0/cuber+a21*fourthr0/fourthr+k21*cuber0/cuber/(1+k22*(1-r0/r)/(1+k23*(1-r0/r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;


    mn=np.zeros((4,4))
    mn[0][0] = -(N2 - W*W*sqsinth) / K2;
    mn[0][3] = -W*r*sqsinth;
    mn[1][1] = Sigma*B*B / N2;
    mn[2][2] = Sigma*sqr;
    mn[3][0] = mn[0][3];
    mn[3][3] = K2*sqr*sqsinth;

    return mn

def metric_KRZ_rderivatives_num(  spin,   d1,   z1,   z2):


    r = z1;
    theta = z2;
    dr = 0.001*r;


    mnm=metric_KRZ(spin, d1, r - dr, theta);
    mnp=metric_KRZ(spin, d1, r + dr, theta);
    
    rdmn=np.zeros((4,4))
    rdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dr;
    rdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dr;
    rdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dr;
    rdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dr;
    rdmn[3][0] = rdmn[0][3];
    rdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dr;

    return rdmn
 

def metric_KRZ_thderivatives_num(  spin,   d1,   z1,   z2):


    r = z1;
    theta = z2;
    dtheta = 0.01;


    metric_KRZ(spin, d1, r, theta - dtheta, mnm);
    metric_KRZ(spin, d1, r, theta + dtheta, mnp);
    
    thdmn=np.zeros((4,4))
    thdmn[0][0] = (mnp[0][0] - mnm[0][0])*0.5 / dtheta;
    thdmn[0][3] = (mnp[0][3] - mnm[0][3])*0.5 / dtheta;
    thdmn[1][1] = (mnp[1][1] - mnm[1][1])*0.5 / dtheta;
    thdmn[2][2] = (mnp[2][2] - mnm[2][2])*0.5 / dtheta;
    thdmn[3][0] = thdmn[0][3];
    thdmn[3][3] = (mnp[3][3] - mnm[3][3])*0.5 / dtheta;

    return thdmn
 

def metric_KRZ_rderivatives(  spin,   d1,   z1,   z2) :

    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = cos(th);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;
    k23 = sqspin / sqr0;

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
    #2018-4-4
    #the oringinal N2 was wrong, see the correct KRZ metric

    rderN2 = (1 - r0 / r)*((2 - r0) / sqr - (2 * (-e0 + k00)*sqr0) / cuber - (3 * d1*cuber0) / fourthr) +    (r0*(1 - (2 - r0) / r + ((-e0 + k00)*sqr0) / sqr + (d1*cuber0) / cuber)) / sqr    +sqcosth*(-3*a20*cuber0/fourthr-4*a21*fourthr0/fourthr/r-3*k21*cuber0/fourthr/ (   1    +      k22*(1-r0/r)/(1+k23*(1-r0/r))  )-     k21*cuber0/cuber/  (1+  k22*(1-r0/r)/(1+k23*(1-r0/r ))  ) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))  )  *  ( k22*r0/sqr/(1+  k23*(1-r0/r)  ) - k22*(1-r0/r) / (1 + k23* (1 - r0 / r)) / (1 + k23* (1 - r0 / r)) *k23*r0/sqr  )    )
    #2018-4-4 modified with the correct metric expression

    '''+
    ((-3 * (-0 + fourthspin / fourthr0)*cuber0) / fourthr - (4 * a21*fourthr0) / np.power(r, 5))*sqcosth''';
    rderB = (-2 * 0 * sqr0) / cuber - (2 * 0 * sqr0*sqcosth) / cuber;
    rderSigma = (-2 * sqcosth*sqspin) / cuber;
    rderW = (2 * sqcosth*((0 * np.power(r0, 3)) / np.power(r, 3) + (2 * spin) / np.power(r, 2) + (0 * np.power(r0, 3)*sqcosth) / np.power(r, 3))*sqspin) /    (np.power(r, 3)*np.power(1 + (sqcosth*sqspin) / np.power(r, 2), 2)) +    ((-3 * 0 * np.power(r0, 3)) / np.power(r, 4) - (4 * spin) / np.power(r, 3) - (3 * 0 * np.power(r0, 3)*sqcosth) / np.power(r, 4)) /    (1 + (sqcosth*sqspin) / np.power(r, 2));
    '''	  rderK2 = (2 * cubespin*sqcosth*((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
    (fourthr*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (2 * sqcosth*sqspin*((k21*cuber0*sqcosth) / cuber + sqspin / sqr)) /
    (cuber*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (spin*((-3 *  0*cuber0) / fourthr - (4 * spin) / cuber - (3 * 0*cuber0*sqcosth) / fourthr)) /
    (r*(1 + (sqcosth*sqspin) / sqr)) - (spin*
    ((0*cuber0) / cuber + (2 * spin) / sqr + (0*cuber0*sqcosth) / cuber)) /
    (sqr*(1 + (sqcosth*sqspin) / sqr)) +
    ((-3 * k21*cuber0*sqcosth) / fourthr - (2 * sqspin) / cuber) / (1 + (sqcosth*sqspin) / sqr);'''
    rderK2 = spin / r*rderW - spin*W / sqr - rderSigma / Sigma / Sigma * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))    + (-2 * k00 * sqr0 / cuber - 3 * k21 * cuber0 * sqcosth / fourthr / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))    - k21 * sqcosth * cuber0 / cuber * (r0 * k22 / sqr / (1 + k23*(1 - r0 / r)) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))) / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26

    rdmn=np.zeros((4,4))
    rdmn[0][0] = (-rderN2 + 2 * rderW*sqsinth*W) / K2 - (rderK2*(-N2 + sqsinth*np.power(W, 2))) / np.power(K2, 2);

    rdmn[1][1] = (np.power(B, 2)*rderSigma) / N2 + (2 * B*rderB*Sigma) / N2 - (np.power(B, 2)*rderN2*Sigma) / np.power(N2, 2);

    rdmn[2][2] = np.power(r, 2)*rderSigma + 2 * r*Sigma;

    rdmn[3][3] = 2 * K2*r*sqsinth + np.power(r, 2)*rderK2*sqsinth;

    rdmn[0][3] = -(r*rderW*sqsinth) - sqsinth*W;

    rdmn[3][0] = rdmn[0][3];

    return rdmn


def metric_KRZ_thderivatives(  spin,   d1,   z1,   z2) :
    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = cos(th); #2018-8-19 之前算cosine居然是开根号。。。不知道以前的raytracing是不是也出过这种问题
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;#2017-10-26
    k23 = sqspin / sqr0;#2017-10-26

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;
    #2018-4-4

    thderN2 = -2 * costh*sinth* (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22* (1 - r0 / r) / (1 + k23 * (1 - r0 / r))));
    #2018-4-4 modified with the correct metric expression
    ''' -2 * costh*(((-0 + fourthspin / fourthr0)*cuber0) / cuber + (a21*fourthr0) / fourthr)*sinth ''';
    thderB = (-2 * costh * 0 * sqr0*sinth) / sqr;
    thderSigma = (-2 * costh*sinth*sqspin) / sqr;
    thderW = (2 * costh*sinth*((0 * cuber0) / cuber + (sqcosth * 0 * cuber0) / cuber + (2 * spin) / sqr)*    sqspin) / (sqr*np.power(1+(sqcosth*sqspin)/sqr, 2)) -    (2 * costh * 0 * cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr));
    '''  thderK2 = (2 * costh*cubespin*sinth*((0*cuber0) / cuber + (sqcosth*0*cuber0) / cuber +
    (2 * spin) / sqr)) / (cuber*np.power(1 + (sqcosth*sqspin) / sqr, 2)) +
    (2 * costh*sinth*sqspin*((sqcosth*k21*cuber0) / cuber + sqspin / sqr)) /
    (sqr*np.power(1 + (sqcosth*sqspin) / sqr, 2)) -
    (2 * costh*k21*cuber0*sinth) / (cuber*(1 + (sqcosth*sqspin) / sqr)) -
    (2 * costh*0*cuber0*sinth*spin) / (fourthr*(1 + (sqcosth*sqspin) / sqr));'''

    thderK2 = spin / r * thderW - thderSigma / Sigma / Sigma  * (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))    + (-2* costh * sinth *  k21 * cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26

    thdmn=np.zeros((4,4))
    
    thdmn[0][0] = (-thderN2 + 2 * np.power(sinth, 2)*thderW*W + 2 * costh*sinth*np.power(W, 2)) / K2 -    (thderK2*(-N2 + np.power(sinth, 2)*np.power(W, 2))) / np.power(K2, 2);
    thdmn[1][1] = (2 * B*Sigma*thderB) / N2 - (np.power(B, 2)*Sigma*thderN2) / np.power(N2, 2) + (np.power(B, 2)*thderSigma) / N2;
    thdmn[2][2] = np.power(r, 2)*thderSigma;
    thdmn[3][3] = 2 * costh*K2*np.power(r, 2)*sinth + np.power(r, 2)*np.power(sinth, 2)*thderK2;
    thdmn[0][3] = -(r*np.power(sinth, 2)*thderW) - 2 * costh*r*sinth*W;
    thdmn[3][0] = thdmn[0][3];


    return thdmn
def Christoffel_KRZ(spin,d1,w1,w2):
    rDg=metric_KRZ_rderivatives(spin, d1,  w1, w2);
    thDg=metric_KRZ_thderivatives(spin, d1, w1, w2);
    Dg=np.zeros((4,4,4))
    Dg[0][0][1]=rDg[0][0];Dg[0][1][1]=rDg[0][1];Dg[0][2][1]=rDg[0][2];Dg[0][3][1]=rDg[0][3];
    Dg[1][0][1]=rDg[1][0];Dg[1][1][1]=rDg[1][1];Dg[1][2][1]=rDg[1][2];Dg[1][3][1]=rDg[1][3];
    Dg[2][0][1]=rDg[2][0];Dg[2][1][1]=rDg[2][1];Dg[2][2][1]=rDg[2][2];Dg[2][3][1]=rDg[2][3];
    Dg[3][0][1]=rDg[3][0];Dg[3][1][1]=rDg[3][1];Dg[3][2][1]=rDg[3][2];Dg[3][3][1]=rDg[3][3];
    Dg[0][0][2]=thDg[0][0];Dg[0][1][2]=thDg[0][1];Dg[0][2][2]=thDg[0][2];Dg[0][3][2]=thDg[0][3];
    Dg[1][0][2]=thDg[1][0];Dg[1][1][2]=thDg[1][1];Dg[1][2][2]=thDg[1][2];Dg[1][3][2]=thDg[1][3];
    Dg[2][0][2]=thDg[2][0];Dg[2][1][2]=thDg[2][1];Dg[2][2][2]=thDg[2][2];Dg[2][3][2]=thDg[2][3];
    Dg[3][0][2]=thDg[3][0];Dg[3][1][2]=thDg[3][1];Dg[3][2][2]=thDg[3][2];Dg[3][3][2]=thDg[3][3];
    g=metric_KRZ(spin, d1, w1, w2)
    invg=metric_KRZ_inverse(spin, d1, w1, w2)
    CS=np.zeros((4,4,4))
    CS[0][0][0] = 0;
    CS[0][0][1] = (invg[0][0] * Dg[0][0][1] + invg[0][3] * Dg[0][3][1]);
    CS[0][0][2] = (invg[0][0] * Dg[0][0][2] + invg[0][3] * Dg[0][3][2]);
    CS[0][0][3] = 0;
    CS[0][1][0] = CS[0][0][1];
    CS[0][1][1] = 0;
    CS[0][1][2] = 0;
    CS[0][1][3] = (invg[0][0] * Dg[0][3][1] + invg[0][3] * Dg[3][3][1]);
    CS[0][2][0] = CS[0][0][2];
    CS[0][2][1] = 0;
    CS[0][2][2] = 0;
    CS[0][2][3] = (invg[0][0] * Dg[0][3][2] + invg[0][3] * Dg[3][3][2]);
    CS[0][3][0] = 0;
    CS[0][3][1] = CS[0][1][3];
    CS[0][3][2] = CS[0][2][3];
    CS[0][3][3] = 0;

    CS[1][0][0] = -invg[1][1] * Dg[0][0][1];
    CS[1][0][1] = 0;
    CS[1][0][2] = 0;
    CS[1][0][3] = -invg[1][1] * Dg[0][3][1];
    CS[1][1][0] = 0;
    CS[1][1][1] = invg[1][1] * Dg[1][1][1];
    CS[1][1][2] = invg[1][1] * Dg[1][1][2];
    CS[1][1][3] = 0;
    CS[1][2][0] = 0;
    CS[1][2][1] = CS[1][1][2];
    CS[1][2][2] = -invg[1][1] * Dg[2][2][1];
    CS[1][2][3] = 0;
    CS[1][3][0] = CS[1][0][3];
    CS[1][3][1] = 0;
    CS[1][3][2] = 0;
    CS[1][3][3] = -invg[1][1] * Dg[3][3][1];

    CS[2][0][0] = -invg[2][2] * Dg[0][0][2];
    CS[2][0][1] = 0;
    CS[2][0][2] = 0;
    CS[2][0][3] = -invg[2][2] * Dg[0][3][2];
    CS[2][1][0] = 0;
    CS[2][1][1] = -invg[2][2] * Dg[1][1][2];
    CS[2][1][2] = invg[2][2] * Dg[2][2][1];
    CS[2][1][3] = 0;
    CS[2][2][0] = 0;
    CS[2][2][1] = CS[2][1][2];
    CS[2][2][2] = invg[2][2] * Dg[2][2][2];
    CS[2][2][3] = 0;
    CS[2][3][0] = CS[2][0][3];
    CS[2][3][1] = 0;
    CS[2][3][2] = 0;
    CS[2][3][3] = -invg[2][2] * Dg[3][3][2];

    CS[3][0][0] = 0;
    CS[3][0][1] = (invg[3][3] * Dg[0][3][1] + invg[3][0] * Dg[0][0][1]);
    CS[3][0][2] = (invg[3][3] * Dg[0][3][2] + invg[3][0] * Dg[0][0][2]);
    CS[3][0][3] = 0;
    CS[3][1][0] = CS[3][0][1];
    CS[3][1][1] = 0;
    CS[3][1][2] = 0;
    CS[3][1][3] = (invg[3][3] * Dg[3][3][1] + invg[3][0] * Dg[0][3][1]);
    CS[3][2][0] = CS[3][0][2];
    CS[3][2][1] = 0;
    CS[3][2][2] = 0;
    CS[3][2][3] = (invg[3][3] * Dg[3][3][2] + invg[3][0] * Dg[0][3][2]);
    CS[3][3][0] = 0;
    CS[3][3][1] = CS[3][1][3];
    CS[3][3][2] = CS[3][2][3];
    CS[3][3][3] = 0;
    return CS

def metric_KRZ_inverse(  spin,   d1,   z1,   z2) :
    r = z1;
    th = z2;

    sqr = r*r;
    cuber = r*sqr;
    fourthr = sqr*sqr;

    sinth = sin(th);
    costh = sqrt(1 - sinth*sinth);
    sqsinth = sinth*sinth;
    sqcosth = costh*costh;

    sqspin = spin*spin;
    cubespin = sqspin*spin;
    fourthspin = sqspin*sqspin;

    r0 = 1 + sqrt(1 - sqspin);
    sqr0 = r0*r0;
    cuber0 = sqr0*r0;
    fourthr0 = sqr0*sqr0;

    a20 = (2 * sqspin) / cuber0;
    a21 = -fourthspin / fourthr0 + 0;
    e0 = (2 - r0) / r0;
    k00 = sqspin / sqr0;
    k21 = fourthspin / fourthr0 - 2 * sqspin / cuber0 - 0;
    w00 = 2 * spin / sqr0;

    k22 = -sqspin / sqr0;#2017-10-26
    k23 = sqspin / sqr0;#2017-10-26
    '''
    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + ((k21 + a20)*cuber0 / cuber + a21*fourthr0 / fourthr)*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26
    '''

    N2 = (1 - r0 / r)*(1 - e0*r0 / r + (k00 - e0)*sqr0 / sqr + d1*cuber0 / cuber) + (a20*cuber0 / cuber + a21*fourthr0 / fourthr + k21*cuber0 / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r))))*sqcosth;
    B = 1 + 0 * sqr0 / sqr + 0 * sqr0*sqcosth / sqr;
    Sigma = 1 + sqspin*sqcosth / sqr;
    W = (w00*sqr0 / sqr + 0 * cuber0 / cuber + 0 * cuber0*sqcosth / cuber) / Sigma;
    K2 = 1 + spin*W / r + (k00*sqr0 / sqr + k21*cuber0*sqcosth / cuber / (1 + k22*(1 - r0 / r) / (1 + k23*(1 - r0 / r)))) / Sigma;#2017-10-26


    invg=np.zeros((4,4))
    invg[0][0] = -(K2 / N2);
    invg[1][1] = N2 / (np.power(B, 2)*Sigma);
    invg[2][2] = 1 / (sqr*Sigma);
    invg[3][3] = (-np.power(W, 2) + N2 / np.power(sin(th), 2)) / (K2*N2*sqr);
    invg[0][3] = -(W / (N2*r));
    invg[3][0] = invg[0][3];


    '''  g[4][4], gg;
    metric_KRZ(spin, d1, z1, z2, g);
    #metric_KRZ_inverse(spin, d1,  w1, w2, invg);
    gg = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    if (std::fabs(gg) > 1e10) :
    invg[0][0] = g[3][3] / gg;
    invg[0][3] = -g[0][3] / gg;
    invg[1][1] = 1 / g[1][1];
    invg[2][2] = 1 / g[2][2];
    invg[3][0] = invg[0][3];
    invg[3][3] = g[0][0] / gg;
    '''
    return invg



def getwave(filename,THETA=np.pi/4,PHI=0,M=1e6,R_pc=5e9,mu=1e-5,usenp=True):
#读入文件名，和观测角，输出引力波
    try:
        index, tau,t,r,th,phi,ut,ur,uth,uphi,F_t,F_r,F_th,F_phi=np.loadtxt(filename,unpack=True)
    except:
        try:
            index, myt_sec, tau,t,r,th,phi,ut,ur,uth,uphi,F_t,F_r,F_th,F_phi=np.loadtxt(filename,unpack=True)
        except:
            print(filename+'  does not exist')
            quit()
    

    if usenp:
        
        x=r*np.sin(th)*np.cos(phi);
        y=r*np.sin(th)*np.sin(phi);
        z=(r*np.cos(th));
        t_tau_dot=(ut)
        x_tau_dot=(ur*np.sin(th)*np.cos(phi) + r*np.cos(th)*np.cos(phi)*uth - r*np.sin(th)*np.sin(phi)*uphi )
        y_tau_dot=(ur*np.sin(th)*np.sin(phi) + r*np.cos(th)*np.sin(phi)*uth + r*np.sin(th)*np.cos(phi)*uphi )
        z_tau_dot=(ur*np.cos(th) - r*np.sin(th)*uth)
        x_t_dot=(x_tau_dot/t_tau_dot)
        y_t_dot=(y_tau_dot/t_tau_dot)
        z_t_dot=(z_tau_dot/t_tau_dot)

        vr_tau_dot=( (F_r*t_tau_dot-ur*F_t)/t_tau_dot/t_tau_dot )
        vth_tau_dot=( (F_th*t_tau_dot-uth*F_t)/t_tau_dot/t_tau_dot )
        vphi_tau_dot=( (F_phi*t_tau_dot-uphi*F_t)/t_tau_dot/t_tau_dot )

        vx_tau_dot=( vr_tau_dot*np.sin(th)*np.cos(phi) + ur/ut*np.cos(th)*np.cos(phi)*uth - ur/ut*np.sin(th)*np.sin(phi)*uphi\
                          + ur*cos(th)*cos(phi)*uth/ut - r*sin(th)*cos(phi)*uth/ut*uth -r*cos(th)*sin(phi)*uth/ut*uphi +r*cos(th)*cos(phi)*vth_tau_dot  \
                          - ur*sin(th)*sin(phi)*uphi/ut - r*cos(th)*sin(phi)*uphi/ut*uth - r*sin(th)*cos(phi)*uphi/ut*uphi - r*sin(th)*sin(phi)*vphi_tau_dot)

        vy_tau_dot=( vr_tau_dot*np.sin(th)*np.sin(phi) + ur/ut*np.cos(th)*np.sin(phi)*uth + ur/ut*np.sin(th)*np.cos(phi)*uphi\
                          + ur*cos(th)*sin(phi)*uth/ut - r*sin(th)*sin(phi)*uth/ut*uth +r*cos(th)*cos(phi)*uth/ut*uphi +r*cos(th)*sin(phi)*vth_tau_dot  \
                          + ur*sin(th)*cos(phi)*uphi/ut + r*cos(th)*cos(phi)*uphi/ut*uth - r*sin(th)*sin(phi)*uphi/ut*uphi + r*sin(th)*cos(phi)*vphi_tau_dot)

        vz_tau_dot=( vr_tau_dot*cos(th) -ur/ut*sin(th)*uth \
                          -ur*sin(th)*uth/ut -r*cos(th)*uth/ut*uth - r*sin(th)*vth_tau_dot )

        x_t_2dot=(vx_tau_dot/ut)
        y_t_2dot=(vy_tau_dot/ut)
        z_t_2dot=(vz_tau_dot/ut)

        #四极矩算法，在trace-reversed gauge的metric


        hbar_xx=(4*(x_t_dot*x_t_dot+x*x_t_2dot))
        hbar_yy=(4*(y_t_dot*y_t_dot+y*y_t_2dot))
        hbar_zz=(4*(z_t_dot*z_t_dot+z*z_t_2dot))
        hbar_xy=(2*(y*x_t_2dot+y_t_2dot*x+2*y_t_dot*x_t_dot))
        hbar_yz=(2*(y*z_t_2dot+y_t_2dot*z+2*y_t_dot*z_t_dot))
        hbar_xz=(2*(z*x_t_2dot+z_t_2dot*x+2*z_t_dot*x_t_dot))

        #由trace-reversed gauge转换到transverse traceless gauge



        hTT_TT=( np.cos(THETA)*np.cos(THETA)* (hbar_xx*np.cos(PHI)*np.cos(PHI) + hbar_xy*np.sin(2*PHI) + hbar_yy*np.sin(PHI)*np.sin(PHI) )  +  hbar_zz*np.sin(THETA)*np.sin(THETA)  -  np.sin(2*THETA)* (hbar_xz*np.cos(PHI)+hbar_yz*np.sin(PHI))  )
        hTT_TP=( np.cos(THETA)* (-0.5*hbar_xx*np.sin(2*PHI) + hbar_xy*np.cos(2*PHI) + 0.5*hbar_yy*np.sin(2*PHI))  +  np.sin(THETA)* (hbar_xz*np.sin(PHI)-hbar_yz*np.cos(PHI)) )
        hTT_PP=( hbar_xx*np.sin(PHI)*np.sin(PHI)  -  hbar_xy*np.sin(2*PHI)  +  hbar_yy*np.cos(PHI)*np.cos(PHI) )
        hTT_plus=(0.5*(hTT_TT-hTT_PP))
        hTT_cross=(hTT_TP)

    else:#不用np用append。。。
        #qseudo_flat spacetime
        x=[];
        y=[];
        z=[];
        t_tau_dot=[]
        z_tau_dot=[]
        y_tau_dot=[]
        x_tau_dot=[]
        z_t_dot=[]
        y_t_dot=[]
        x_t_dot=[]
        vr_tau_dot=[]
        vth_tau_dot=[]
        vphi_tau_dot=[]
        vx_tau_dot=[]
        vy_tau_dot=[]
        vz_tau_dot=[]
        x_t_2dot=[]
        y_t_2dot=[]
        z_t_2dot=[]
        for i in np.arange(index.size):
            x.append(r[i]*np.sin(th[i])*np.cos(phi[i]));
            y.append(r[i]*np.sin(th[i])*np.sin(phi[i]));
            z.append(r[i]*np.cos(th[i]));
            t_tau_dot.append(ut[i])
            x_tau_dot.append(ur[i]*np.sin(th[i])*np.cos(phi[i]) + r[i]*np.cos(th[i])*np.cos(phi[i])*uth[i] - r[i]*np.sin(th[i])*np.sin(phi[i])*uphi[i] )
            y_tau_dot.append(ur[i]*np.sin(th[i])*np.sin(phi[i]) + r[i]*np.cos(th[i])*np.sin(phi[i])*uth[i] + r[i]*np.sin(th[i])*np.cos(phi[i])*uphi[i] )
            z_tau_dot.append(ur[i]*np.cos(th[i]) - r[i]*np.sin(th[i])*uth[i])
            x_t_dot.append(x_tau_dot[i]/t_tau_dot[i])
            y_t_dot.append(y_tau_dot[i]/t_tau_dot[i])
            z_t_dot.append(z_tau_dot[i]/t_tau_dot[i])

            vr_tau_dot.append( (F_r[i]*t_tau_dot[i]-ur[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )
            vth_tau_dot.append( (F_th[i]*t_tau_dot[i]-uth[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )
            vphi_tau_dot.append( (F_phi[i]*t_tau_dot[i]-uphi[i]*F_t[i])/t_tau_dot[i]/t_tau_dot[i] )

            vx_tau_dot.append( vr_tau_dot[i]*np.sin(th[i])*np.cos(phi[i]) + ur[i]/ut[i]*np.cos(th[i])*np.cos(phi[i])*uth[i] - ur[i]/ut[i]*np.sin(th[i])*np.sin(phi[i])*uphi[i]\
                 + ur[i]*cos(th[i])*cos(phi[i])*uth[i]/ut[i] - r[i]*sin(th[i])*cos(phi[i])*uth[i]/ut[i]*uth[i] -r[i]*cos(th[i])*sin(phi[i])*uth[i]/ut[i]*uphi[i] +r[i]*cos(th[i])*cos(phi[i])*vth_tau_dot[i]  \
                 - ur[i]*sin(th[i])*sin(phi[i])*uphi[i]/ut[i] - r[i]*cos(th[i])*sin(phi[i])*uphi[i]/ut[i]*uth[i] - r[i]*sin(th[i])*cos(phi[i])*uphi[i]/ut[i]*uphi[i] - r[i]*sin(th[i])*sin(phi[i])*vphi_tau_dot[i])

            vy_tau_dot.append( vr_tau_dot[i]*np.sin(th[i])*np.sin(phi[i]) + ur[i]/ut[i]*np.cos(th[i])*np.sin(phi[i])*uth[i] + ur[i]/ut[i]*np.sin(th[i])*np.cos(phi[i])*uphi[i]\
                 + ur[i]*cos(th[i])*sin(phi[i])*uth[i]/ut[i] - r[i]*sin(th[i])*sin(phi[i])*uth[i]/ut[i]*uth[i] +r[i]*cos(th[i])*cos(phi[i])*uth[i]/ut[i]*uphi[i] +r[i]*cos(th[i])*sin(phi[i])*vth_tau_dot[i]  \
                 + ur[i]*sin(th[i])*cos(phi[i])*uphi[i]/ut[i] + r[i]*cos(th[i])*cos(phi[i])*uphi[i]/ut[i]*uth[i] - r[i]*sin(th[i])*sin(phi[i])*uphi[i]/ut[i]*uphi[i] + r[i]*sin(th[i])*cos(phi[i])*vphi_tau_dot[i])

            vz_tau_dot.append( vr_tau_dot[i]*cos(th[i]) -ur[i]/ut[i]*sin(th[i])*uth[i] \
                             -ur[i]*sin(th[i])*uth[i]/ut[i] -r[i]*cos(th[i])*uth[i]/ut[i]*uth[i] - r[i]*sin(th[i])*vth_tau_dot[i] )

            x_t_2dot.append(vx_tau_dot[i]/ut[i])
            y_t_2dot.append(vy_tau_dot[i]/ut[i])
            z_t_2dot.append(vz_tau_dot[i]/ut[i])

        #四极矩算法，在trace-reversed gauge的metric

        hbar_xx=[]
        hbar_yy=[]
        hbar_zz=[]
        hbar_xy=[]
        hbar_yz=[]
        hbar_xz=[]
        for i in np.arange(index.size):
            hbar_xx.append(4*(x_t_dot[i]*x_t_dot[i]+x[i]*x_t_2dot[i]))
            hbar_yy.append(4*(y_t_dot[i]*y_t_dot[i]+y[i]*y_t_2dot[i]))
            hbar_zz.append(4*(z_t_dot[i]*z_t_dot[i]+z[i]*z_t_2dot[i]))
            hbar_xy.append(2*(y[i]*x_t_2dot[i]+y_t_2dot[i]*x[i]+2*y_t_dot[i]*x_t_dot[i]))
            hbar_yz.append(2*(y[i]*z_t_2dot[i]+y_t_2dot[i]*z[i]+2*y_t_dot[i]*z_t_dot[i]))
            hbar_xz.append(2*(z[i]*x_t_2dot[i]+z_t_2dot[i]*x[i]+2*z_t_dot[i]*x_t_dot[i]))

        #由trace-reversed gauge转换到transverse traceless gauge

        hTT_TT=[]
        hTT_PP=[]
        hTT_TP=[]
        hTT_plus=[]
        hTT_cross=[]

        for i in np.arange(index.size):


            hTT_TT.append( np.cos(THETA)*np.cos(THETA)* (hbar_xx[i]*np.cos(PHI)*np.cos(PHI) + hbar_xy[i]*np.sin(2*PHI) + hbar_yy[i]*np.sin(PHI)*np.sin(PHI) )  +  hbar_zz[i]*np.sin(THETA)*np.sin(THETA)  -  np.sin(2*THETA)* (hbar_xz[i]*np.cos(PHI)+hbar_yz[i]*np.sin(PHI))  )
            hTT_TP.append( np.cos(THETA)* (-0.5*hbar_xx[i]*np.sin(2*PHI) + hbar_xy[i]*np.cos(2*PHI) + 0.5*hbar_yy[i]*np.sin(2*PHI))  +  np.sin(THETA)* (hbar_xz[i]*np.sin(PHI)-hbar_yz[i]*np.cos(PHI)) )
            hTT_PP.append( hbar_xx[i]*np.sin(PHI)*np.sin(PHI)  -  hbar_xy[i]*np.sin(2*PHI)  +  hbar_yy[i]*np.cos(PHI)*np.cos(PHI) )
            hTT_plus.append(0.5*(hTT_TT[i]-hTT_PP[i]))
            hTT_cross.append(hTT_TP[i])

    #注意上面算出来的h还要*mu（mass ratio）/R（观测距离，也以M为单位）才是真的strain
    #发现一个小问题，上面定义的数据类型大部分都是list，但是array才比较好用
    #还要注意一点几何单位制和SI单位的转换

    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #M=1e6 # clight*clight*clight/Grav/Msol/1 #中心天体质量，以太阳质量为单位

    #把时间转换成秒
    t_sec=t*M*Msol*Grav/clight/clight/clight
    #dt=t_sec[1]-t_sec[0]

    #把pc距离转换成M为单位
    #R_pc=5e9  #以pc为单位的观测距离
    R=R_pc*3.0857e16*clight*clight/Grav/M/Msol  #以中心天体质量为单位的，长度米与中心天体质量的换算是 1m/kg = clight*clight/G

    #小天体的质量
    #mu=1e-5 #应该是以中心天体质量为单位的

    hTT_plus_true=np.array(hTT_plus)*mu/R*M
    hTT_cross_true=np.array(hTT_cross)*mu/R*M

    ########用于计算的波形，plus作为实部，cross作为虚部

    return t_sec,hTT_plus_true+hTT_cross_true*1j
    
def bracket(mydata,mytemp,dt,fnoise=[-1,-1],Snoise=[-1,-1]):
#算mydata和mytemp的内积
#fnoise-Snoise是噪声的功率谱
    
    #两个时间序列长度要一样，不一样就取最小的
    if len(mydata) != len(mytemp):
        print("inner product: length not match!")

    #采样频率
    fs=1/dt

    #做傅里叶变换，注意要除以采样频率才是真的amplitude
    mydata_fft=np.fft.fft(mydata)/len(mydata)
    mytemp_fft=np.fft.fft(mytemp)/len(mydata)
    #变换后的频率序列,注意，只有前一半是正的频率
    freq=np.fft.fftfreq(min( len(mydata),len(mytemp)),dt)
    
    #默认为LISA noise
    if fnoise[0]==-1:
        ##########LISA noise, reference: https://arxiv.org/abs/gr-qc/0607007v2
        u=2*np.pi*freq*50/3 #见reference（36）上面一段
        Sn=[]  #LISA noise
        for i in np.arange(freq.size/2):
            i=int(i)
            if i==0:
                Sn.append(1e10)
            elif(u[i]<0.25):
                Sn.append(8.08e-48/((2*np.pi*freq[i])**4) +5.52e-41 )
            else :
                Sn.append( (2.88e-48/((2*np.pi*freq[i])**4) +5.52e-41 ) *u[i]*u[i]/ ( (1+cos(u[i])*cos(u[i]) )*(1.0/3.0-2.0/u[i]/u[i]) + sin(u[i])**2 + 4*sin(u[i])*cos(u[i])/(u[i]**3) ) )

    else:
        try:#如果输入了噪声功率谱,先直接interpolate
            noisefit=interp1d(fnoise,Snoise)
            Sn=noisefit(freq[0:int(freq.size/2)])#如果输入的范围不够在这行会报错
            
        except:#如果报错说明输入的功率谱范围覆盖不了信号傅里叶变换后的功率谱，那么把输入功率谱的左右两侧用直线外推
            leftfit=np.poly1d(np.polyfit(fnoise[0:int(freq.size/8)],Snoise[0:int(freq.size/8)],1))
            rightfit=np.poly1d(np.polyfit(fnoise[int(freq.size/2)-int(freq.size/8):int(freq.size/2)],Snoise[int(freq.size/2)-int(freq.size/8):int(freq.size/2)],1))
            Sn=[]
            for i in np.arange(int(freq.size/2)):
                if freq[i]<fnoise[0]:
                    Sn.append(leftfit(freq[i]))
                elif freq[i]>fnoise[-1]:
                    Sn.append(rightfit(freq[i]))
                else:
                    Sn.append(noisefit(freq[i]))
                    
    ###########SNR
    SNR=0
    df=freq[1]-freq[0]
    for i in np.arange(int(freq.size/2)):
        i=int(i)
        SNR=SNR+(( mydata_fft[i]*np.conjugate(mytemp_fft[i]) + np.conjugate(mydata_fft[i])*mytemp_fft[i])/(Sn[i]))*df
    
    return np.abs(SNR)

def getfreq_fromtraj(tau,r,phi):
    #由序列获得orbital frequency
    p=np.mean(r)
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(tau.size-1):
        if(r[i]>p and r[i+1]<=p):
            indr.append(i)

    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(tau[indr[ii]]-tau[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(tau[indr[ii]]-tau[indr[ii-1]]))

    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))

    return avgomgr,avgomgphi
def getfreq_frommaxi(tau,r,phi):
    #由序列获得orbital frequency,这个序列必须是从r最大值开始的
    #p=np.mean(r)
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(tau.size-1):
        if i==0:
            indr.append(i)
        elif(r[i]>r[i-1] and r[i+1]<r[i]):
            indr.append(i)

    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(tau[indr[ii]]-tau[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(tau[indr[ii]]-tau[indr[ii-1]]))

    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))

    return avgomgr,avgomgphi

def getfreq_dt_fromepa_npsum(e,p,spin):
    if np.abs(e-0.0)>1e-8:
        rmax=p/(1-e)
        rmin=p/(1+e)
        invgmin=metric_KRZ_inverse(spin,0,rmin,np.pi/2)
        invgmax=metric_KRZ_inverse(spin,0,rmax,np.pi/2)

        EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / ( invgmax[0][0]-invgmin[0][0] );
        Lz = sqrt( (invgmax[3][0]-invgmin[3][0]) / ( EoverL*EoverL*( invgmin[3][0]*invgmax[0][0] - invgmax[3][0]*invgmin[0][0] )+ ( invgmin[3][0]*invgmax[3][3]- invgmax[3][0]*invgmin[3][3] )  )   );
        E=Lz*EoverL
        
    else:
        Gamma=Christoffel_KRZ(spin, 0, p, th);
        g=metric_KRZ(spin, 0, p, th);
        utoverup = (-Gamma[1][0][3] + np.sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
        testup = np.sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
        testut = testup*utoverup;
        E = -g[0][0] * testut - g[0][3] * testup;
        Lz = g[0][3] * testut + g[3][3] * testup;
        
    x=Lz-spin*E
    
    dchi=1e-6
    chi=np.linspace(dchi/2,np.pi-dchi/2,int(np.pi/dchi))
    J=1-2*(1+e*np.cos(chi))/p+spin**2/p**2*(1+e*np.cos(chi))**2
    Vr=x**2+spin**2+2*spin*x*E-2*x**2/p*(3+e*np.cos(chi))
    Vphi=x+spin*E-2*x/p*(1+e*cos(chi))
    Vt=spin**2*E-2*spin*x/p*(1+e*np.cos(chi)) + E*p**2/(1+e*np.cos(chi))**2
    
    Tr=2*np.sum(Vt/J/np.sqrt(Vr))*dchi
    Dphi=2*np.sum(Vphi/J/np.sqrt(Vr))*dchi
    omg_rdt=2*np.pi/Tr
    omg_phidt=Dphi/Tr
    
    return omg_rdt,omg_phidt

def Tr_int(chi,spin,e,p):
    if np.abs(e-0.0)>1e-8:
        rmax=p/(1-e)
        rmin=p/(1+e)
        invgmin=metric_KRZ_inverse(spin,0,rmin,np.pi/2)
        invgmax=metric_KRZ_inverse(spin,0,rmax,np.pi/2)

        EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / ( invgmax[0][0]-invgmin[0][0] );
        Lz = sqrt( (invgmax[3][0]-invgmin[3][0]) / ( EoverL*EoverL*( invgmin[3][0]*invgmax[0][0] - invgmax[3][0]*invgmin[0][0] )+ ( invgmin[3][0]*invgmax[3][3]- invgmax[3][0]*invgmin[3][3] )  )   );
        E=Lz*EoverL
        
    else:
        Gamma=Christoffel_KRZ(spin, 0, p, np.pi/2);
        g=metric_KRZ(spin, 0, p, np.pi/2);
        utoverup = (-Gamma[1][0][3] + np.sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
        testup = np.sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
        testut = testup*utoverup;
        E = -g[0][0] * testut - g[0][3] * testup;
        Lz = g[0][3] * testut + g[3][3] * testup;
    
    x=Lz-spin*E
    J=1-2*(1+e*np.cos(chi))/p+spin**2/p**2*(1+e*np.cos(chi))**2
    Vr=x**2+spin**2+2*spin*x*E-2*x**2/p*(3+e*np.cos(chi))
    Vphi=x+spin*E-2*x/p*(1+e*cos(chi))
    Vt=spin**2*E-2*spin*x/p*(1+e*np.cos(chi)) + E*p**2/(1+e*np.cos(chi))**2

    Tr_=2*Vt/J/np.sqrt(Vr)
    return Tr_
def Dphi_int(chi,spin,e,p):
    if np.abs(e-0.0)>1e-8:
        rmax=p/(1-e)
        rmin=p/(1+e)
        invgmin=metric_KRZ_inverse(spin,0,rmin,np.pi/2)
        invgmax=metric_KRZ_inverse(spin,0,rmax,np.pi/2)

        EoverL = ((invgmax[3][0] - invgmin[3][0]) + sqrt((invgmax[3][0] - invgmin[3][0]) *(invgmax[3][0] - invgmin[3][0]) - (invgmax[0][0] - invgmin[0][0])*(invgmax[3][3] - invgmin[3][3]))) / ( invgmax[0][0]-invgmin[0][0] );
        Lz = sqrt( (invgmax[3][0]-invgmin[3][0]) / ( EoverL*EoverL*( invgmin[3][0]*invgmax[0][0] - invgmax[3][0]*invgmin[0][0] )+ ( invgmin[3][0]*invgmax[3][3]- invgmax[3][0]*invgmin[3][3] )  )   );
        E=Lz*EoverL
        
    else:
        Gamma=Christoffel_KRZ(spin, 0, p, np.pi/2);
        g=metric_KRZ(spin, 0, p, np.pi/2);
        utoverup = (-Gamma[1][0][3] + np.sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
        testup = np.sqrt(-1 / (utoverup*utoverup*g[0][0] + 2 * utoverup*g[0][3] + g[3][3]));
        testut = testup*utoverup;
        E = -g[0][0] * testut - g[0][3] * testup;
        Lz = g[0][3] * testut + g[3][3] * testup;
    x=Lz-spin*E
    J=1-2*(1+e*np.cos(chi))/p+spin**2/p**2*(1+e*np.cos(chi))**2
    Vr=x**2+spin**2+2*spin*x*E-2*x**2/p*(3+e*np.cos(chi))
    Vphi=x+spin*E-2*x/p*(1+e*cos(chi))
    Vt=spin**2*E-2*spin*x/p*(1+e*np.cos(chi)) + E*p**2/(1+e*np.cos(chi))**2

    Dp_=2*Vphi/J/np.sqrt(Vr)
    return Dp_

def getfreq_dt_fromepa(e,p,spin):
    myTr,err=scipy.integrate.quad(Tr_int,0,np.pi,args=(spin,e,p))
    myDphi,err=scipy.integrate.quad(Dphi_int,0,np.pi,args=(spin,e,p))
    omg_rdt=2*np.pi/myTr
    omg_phidt=myDphi/myTr
    return omg_rdt,omg_phidt

def getfreq_dt_frommaxi(t,r,phi):
    #由序列获得orbital frequency(对t的),这个序列必须是从r最大值开始的
    #p=np.mean(r)
    #2018-10-10 最大值处如果两个数正好在第十位小数都相等，会出问题（这个点会没算进去），所以稍微改了一下
    tol=9e-11 #数据里r的精度是1e-10
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(t.size-1):
        if i==0:
            indr.append(i)
        elif(r[i]-r[i-1]>tol and r[i+1]-r[i]<-tol):
            indr.append(i)
############2018-10-10##############
        elif(r[i]-r[i-1]>tol and np.abs(r[i+1]-r[i])<tol  ):
            testi=[]
            for jj in range(100):
                testi.append(i+jj)
                if r[i+jj+1]-r[i+jj]<-tol:
                    break

            indr.append(int(round(np.mean(np.array(testi) ) ) ) )
############2018-10-10##################
    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(t[indr[ii]]-t[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(t[indr[ii]]-t[indr[ii-1]]))

    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))

    return avgomgr,avgomgphi

def getfreq_sec_fromepma(e,p,M,spin):

    omg=np.array(getfreq_dt_fromepa(e,p,spin))
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgsec=omg*clight**3/M/Msol/Grav
    return omgsec
def getfreq_sec_frommaxi(t,r,phi,M):

    omgavg=np.array(getfreq_dt_frommaxi(t,r,phi))
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgavgsec=omgavg*clight**3/M/Msol/Grav
    return omgavgsec

def circfreq_sec_fromrma(r,M,spin):
    Gamma=Christoffel_KRZ(spin, 0, r, np.pi/2);
    g=metric_KRZ(spin, 0, r, np.pi/2);
    utoverup = (-Gamma[1][0][3] + np.sqrt(Gamma[1][0][3] * Gamma[1][0][3] - Gamma[1][3][3] * Gamma[1][0][0])) / Gamma[1][0][0];
    omg=1/utoverup
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgsec=omg*clight**3/M/Msol/Grav
    return omgsec

def circfreq_sec_fromtrace(t,phi,M):
    
    omg=(phi[-1]-phi[0])/(t[-1]-t[0])
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgsec=omg*clight**3/M/Msol/Grav
    return omgsec

########### below: 3D frequency & subroutines
def Vr_freq3sub(E,L,e,p,iota,spin):
    #Vr at apestron and periastron, assuming Vth(theta_min)=0, i.e.Q=cos(theta)*(a**2*(1-E**2)+L**2/(sin(theta)*sin(theta)));
    a=spin
    ra=p/(1-e);
    rp=p/(1+e);
    theta=np.pi/2-iota;

    Del=ra**2-2*ra+a**2;
    P=E*(ra**2+a**2)-a*L;
    Q=np.cos(theta)**2*(a**2*(1-E**2)+L**2/(np.sin(theta)*np.sin(theta)));
    vra=P**2-Del*(ra**2+(L-a*E)**2+Q);

    Delp=rp**2-2*rp+a**2;
    Pp=E*(rp**2+a**2)-a*L;
    vrp=Pp**2-Delp*(rp**2+(L-a*E)**2+Q);
    return vra,vrp

def Wtilt_int(chi,e,p,E,Lz,Q,spin):
    myF=E+spin**2*E/p**2*(1+e*np.cos(chi))**2 - 2*spin*(Lz-spin*E)*(1+e*np.cos(chi))**3 / p**3 
    myJ=(1-E**2)*(1-e**2)+2*(1-E**2-(1-e**2)/p)*(1+e*np.cos(chi))+((1-E**2)*(3+e**2)/(1-e**2) -4/p +(spin**2*(1-E**2)+Lz**2+Q)*(1-e**2)/p**2 )*(1+e*np.cos(chi))**2
    myH=1-2/p*(1+e*cos(chi))+(spin**2)/(p**2)*(1+e*np.cos(chi))**2
    myW_=p**2*myF/( (1+e*np.cos(chi))**2 * myH * np.sqrt(myJ) )
    return myW_

def Xtilt_int(chi,e,p,E,Lz,Q,spin):
    
    myJ=(1-E**2)*(1-e**2)+2*(1-E**2-(1-e**2)/p)*(1+e*np.cos(chi))+((1-E**2)*(3+e**2)/(1-e**2) -4/p +(spin**2*(1-E**2)+Lz**2+Q)*(1-e**2)/p**2 )*np.power(1+e*np.cos(chi),2)
    myX_=1/np.sqrt(myJ)
    return myX_

def Ytilt_int(chi,e,p,E,Lz,Q,spin):
    
    myJ=(1-E**2)*(1-e**2)+2*(1-E**2-(1-e**2)/p)*(1+e*np.cos(chi))+((1-E**2)*(3+e**2)/(1-e**2) -4/p +(spin**2*(1-E**2)+Lz**2+Q)*(1-e**2)/p**2 )*(1+e*np.cos(chi))**2
    myH=1-2/p*(1+e*cos(chi))+(spin**2)/(p**2)*(1+e*np.cos(chi))**2
    myG=Lz-2*(Lz-spin*E)/p*(1+e*np.cos(chi))

    myY_=p**2 /(1+e*np.cos(chi))**2/ np.sqrt(myJ) 
    return myY_

def Ztilt_int(chi,e,p,E,Lz,Q,spin):
    
    myJ=(1-E**2)*(1-e**2)+2*(1-E**2-(1-e**2)/p)*(1+e*np.cos(chi))+((1-E**2)*(3+e**2)/(1-e**2) -4/p +(spin**2*(1-E**2)+Lz**2+Q)*(1-e**2)/p**2 )*np.power(1+e*np.cos(chi),2)
    myH=1-2/p*(1+e*cos(chi))+(spin**2)/(p**2)*np.power(1+e*np.cos(chi),2)
    myG=Lz-2*(Lz-spin*E)/p*(1+e*np.cos(chi))
    myZ_=np.multiply(myG, np.multiply( np.power(myH,-1),np.power(myJ,-0.5) ) )
    return myZ_

def freq3_dt(e,p,iota,spin):
    def myfunc(x):
        E=x[0]
        L=x[1]
        V= Vr_freq3sub(E,L,e,p,iota,spin)
        #print('E:%f, L:%f, V: %f, %f'%(E,L,V[0],V[1]))
        return V
    E,Lz=fsolve(myfunc,[0.9,3.0])
    th_min=np.pi/2-iota
    Q=np.cos(th_min)**2*(spin**2*(1-E**2)+Lz**2/(np.sin(th_min)*np.sin(th_min)))
    z_plus=np.sqrt(( Q+Lz**2+spin**2*(1-E**2) + np.sqrt( (Q+Lz**2+spin**2*(1-E**2))**2 - 4*Q*spin**2*(1-E**2) ) )/(2*spin**2*(1-E**2)))
    z_minus=np.sqrt(( Q+Lz**2+spin**2*(1-E**2) - np.sqrt( (Q+Lz**2+spin**2*(1-E**2))**2 - 4*Q*spin**2*(1-E**2) ) )/(2*spin**2*(1-E**2)))
    if np.abs(iota-0.0)<1e-20:
        z_minus=0.0
    beta=np.sqrt(spin**2*(1-E**2))
    k=(z_minus/z_plus)**2
    def K_int(psi,k):
        return 1.0/np.sqrt(1-k*np.sin(psi)**2)
    def E_int(psi,k):
        return np.sqrt(1-k*np.sin(psi)**2)
    def PI_int(psi,z_minus,k):
        return 1.0/( (1-z_minus**2*np.sin(psi)**2)*np.sqrt(1-k*np.sin(psi)**2) )
    K,err=scipy.integrate.quad(K_int,0,np.pi/2,args=(k))
    E_k,err=scipy.integrate.quad(E_int,0,np.pi/2,args=(k))
    PI,err=scipy.integrate.quad(PI_int,0,np.pi/2,args=(z_minus,k))
    
    Ytilt,err=scipy.integrate.quad(Ytilt_int,0,np.pi,args=(e,p,E,Lz,Q,spin))
    Ztilt,err=scipy.integrate.quad(Ztilt_int,0,np.pi,args=(e,p,E,Lz,Q,spin))
    Xtilt,err=scipy.integrate.quad(Xtilt_int,0,np.pi,args=(e,p,E,Lz,Q,spin))

    LAMBDA=(Ytilt+spin**2 * z_plus**2*Xtilt)*K-spin**2 *z_plus**2 *Xtilt*E_k

    omgr=np.pi*p*K/(1-e**2)/LAMBDA
    omgth=np.pi*beta*z_plus*Xtilt/2/LAMBDA
    omgphi=( (Ztilt-Lz*Xtilt)*K + Lz*Xtilt*PI )/LAMBDA
    
    Wtilt,err=scipy.integrate.quad(Wtilt_int,0,np.pi,args=(e,p,E,Lz,Q,spin))
    gamma=( (Wtilt+spin**2*z_plus**2*E*Xtilt)*K - spin**2*z_plus**2*E*Xtilt*E_k )/LAMBDA
    
    omgr_dt=omgr/gamma
    omgth_dt=omgth/gamma
    omgphi_dt=omgphi/gamma
    return omgr_dt,omgth_dt,omgphi_dt

def freq3_sec(e,p,iota,spin,M):
    omg=freq3_dt(e,p,iota,spin)
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgsec=np.array(omg)*clight**3/M/Msol/Grav
    return omgsec
######### above: 3D frequency & subroutines
def getELQ(e,p,iota,spin):
    def myfunc(x):
        E=x[0]
        L=x[1]
        V= Vr_freq3sub(E,L,e,p,iota,spin)
        #print('E:%f, L:%f, V: %f, %f'%(E,L,V[0],V[1]))
        return V
    E,Lz=fsolve(myfunc,[0.9,3.0])
    th_min=np.pi/2-iota
    Q=np.cos(th_min)**2*(spin**2*(1-E**2)+Lz**2/(np.sin(th_min)*np.sin(th_min)))
    return E,Lz,Q

def freq3_dt_fromtrace(t,r,th,phi):
    #由序列获得orbital frequency(对t的),这个序列必须是从r最大值开始的
    #p=np.mean(r)
    #2018-10-10 最大值处如果两个数正好在第十位小数都相等，会出问题（这个点会没算进去），所以稍微改了一下
    tol=9e-11 #数据里r的精度是1e-10
    omgr=[]
    omgphi=[]
    indr=[]
    phi=phi-phi[0]
    n=1
    for i in np.arange(t.size-1):
        if i==0:
            indr.append(i)
        elif(r[i]-r[i-1]>tol and r[i+1]-r[i]<-tol):
            indr.append(i)
    ############2018-10-10##############
        elif(r[i]-r[i-1]>tol and np.abs(r[i+1]-r[i])<tol  ):
            testi=[]
            for jj in range(100):
                testi.append(i+jj)
                if r[i+jj+1]-r[i+jj]<-tol:
                    break

            indr.append(int(round(np.mean(np.array(testi) ) ) ) )
    ############2018-10-10##################

    #####2018-10-29 theta的ind
    indth=[]
    omgth=[]
    for i in np.arange(t.size-1):
        if i==0:
            indth.append(i)
        elif th[i]<=np.pi/2 and th[i+1]>np.pi/2:
            indth.append(i)

    for ii in np.arange(len(indr)):
        if ii==0:
            continue
        omgr.append(2*np.pi/(t[indr[ii]]-t[indr[ii-1]]))
        omgphi.append((phi[indr[ii]]-phi[indr[ii-1]])/(t[indr[ii]]-t[indr[ii-1]]))

    #indth=indth[0:799]
    for jj in np.arange(len(indth)):
        if jj==0:
            continue
        omgth.append(2.0*np.pi/(t[indth[jj]]-t[indth[jj-1]]))
    omgr=np.array(omgr)
    avgomgr=np.mean(omgr)
    avgomgphi=np.mean(np.array(omgphi))
    avgomgth=np.mean(np.array(omgth))

    #2018-10-29:试试直接数周期来当做omg
    numofcyc=len(indr)-1
    period=(t[indr[-1]]-t[indr[0]])/float(numofcyc)
    dphi=(phi[indr[-1]]-phi[indr[0]])/float(numofcyc)
    newomgr=2.0*np.pi/period
    newomgphi=dphi/period
    numofcycth=len(indth)-1
    periodth=(t[indth[-1]]-t[indth[0]])/float(numofcycth)
    newomgth=2.0*np.pi/periodth
    
    #2018-11-3:对omgth做个修正，把尾部波动很大不在平均值附近的去掉，再重新数周期
    for i in np.arange(100):
        if i==0:
            continue
        if np.abs(omgth[-i]-newomgth)<2e-3:#精度2e-3,再除以约2e3个周期数，最后omg误差大约只有1e-6，经过1e5秒，还行吧
            break
    indth=indth[0:-i]
    numofcycth=len(indth)-1
    periodth=(t[indth[-1]]-t[indth[0]])/float(numofcycth)
    newomgth=2.0*np.pi/periodth
    return newomgr,newomgth,newomgphi

def freq3_sec_fromtrace(t,r,th,phi,M):

    omg=np.array(freq3_dt_fromtrace(t,r,th,phi))
    ########转换单位
    Grav=6.674e-11 #引力常数
    clight=2.998e8 #光速
    Msol=1.989e30  #太阳质量，以千克做单位

    #把频率换成s^-1
    omgavgsec=omg*clight**3/M/Msol/Grav
    return omgavgsec

def bracket_interp(wave1,wave2,tottime=-1,dt=5.0):
#wave1,wave2满足wave[0]是时间，wave[1]是波形（复数，hplus+i*hcross）,也就是getwave给出的
#LISA的采样间隔是5.0s
    intp1=interp1d(wave1[0],wave1[1],kind='cubic')
    intp2=interp1d(wave2[0],wave2[1],kind='cubic')
    if tottime==-1:#默认是能取多少重叠就取多少重叠
        tottime=min(wave2[0][-1],wave1[0][-1])
    timeseries=np.arange(0,tottime,dt)
    return bracket(intp2(timeseries),intp1(timeseries),dt)

def overlap(wave1,wave2,tottime=-1,dt=5.0):
    if tottime==-1:#默认是能取多少重叠就取多少重叠
        tottime=min(wave2[0][-1],wave1[0][-1])
    return bracket_interp(wave1,wave2,tottime=tottime,dt=dt)/np.sqrt( bracket_interp(wave1,wave1,tottime=tottime,dt=dt) *bracket_interp(wave2,wave2,tottime=tottime,dt=dt) )
