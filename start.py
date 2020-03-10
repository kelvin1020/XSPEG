# coding: utf-8
import numpy
from numpy import sin,cos,sqrt
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import scipy

#z1就是r，z2就是theta

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

    hTT_plus_true=np.array(hTT_plus)*mu/R
    hTT_cross_true=np.array(hTT_cross)*mu/R

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
            leftfit=np.poly1d(np.polyfit(fnoise[0:int(fnoise.size/8)],Snoise[0:int(fnoise.size/8)],1))
            rightfit=np.poly1d(np.polyfit(fnoise[int(fnoise.size)-int(fnoise.size/8):int(fnoise.size)],Snoise[int(fnoise.size)-int(fnoise.size/8):int(fnoise.size)],1))
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
    #reference: PHYSICAL REVIEW D 66, 044002 (2002)
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

###########--------- below: 3D frequency & subroutines
def Vr_freq3sub(E,L,e,p,iota,spin):
    #ref: eq.(30) in PHYSICAL REVIEW D 96, 044005 (2017)
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
#########---------------- above: 3D frequency & subroutines
def RRprime_circ(E,Lz,r,th,a): #anxillary function for getting ELQ for circular orbits
    Delta = r**2 - 2 *r + a**2
    Q = np.cos(th)**2 *(a**2 * (1 - E**2) + Lz**2 /np.sin(th)**2 ) #令Vth=0
    R = (E*(r**2 + a**2) - a*Lz)**2 - Delta*(r**2 + (Lz - a*E)**2 + Q) 
    Rprime= (-(-2 + 2*r))*(((-a)*E + Lz)**2 + Q + r**2) - 2*r*(a**2 - 2*r + r**2) + 4*E*r*((-a)*Lz + E*(a**2 + r**2))
    return R,Rprime

def getELQ(e,p,iota,spin):
    if(e<1e-5):#circular orbits
        r=p
        th=np.pi/2-iota
        def myfunc(x):
            E=x[0]
            Lz=x[1]
            return RRprime_circ(E,Lz,r,th,spin)
        output=fsolve(myfunc,[0.9,3.0],full_output=True,xtol=1e-20)
        E=output[0][0];Lz=output[0][1]
        Q = np.cos(th)**2 *(spin**2 * (1 - E**2) + Lz**2 /np.sin(th)**2 )
    else:
        def myfunc(x):
            E=x[0]
            L=x[1]
            V= Vr_freq3sub(E,L,e,p,iota,spin)
            #print('E:%f, L:%f, V: %f, %f'%(E,L,V[0],V[1]))
            return V
        E,Lz=fsolve(myfunc,[0.9,3.0])
        th_min=np.pi/2-iota
        Q=np.cos(th_min)**2*(spin**2*(1-E**2)+Lz**2/(np.sin(th_min)*np.sin(th_min)))
    return np.array([E,Lz,Q])

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


# XSPEGdir=os.environ['XSPEGLIB']#############################################################################################
XSPEGdir = './local'


class readINCAR:
    def __init__(self, filename):
        for line in open(filename):
            if (line[0:15].strip() != '' and line.strip()[0]!='#'):
                key = line[0:15].strip()
                mon = line[16:].strip()
                setattr(self, key.lower(), eval(mon))


#initialization

curdir=os.environ['PWD'];
try:
    incarfile=open(curdir+'/INCAR')
    outcarfile=open(curdir+'/OUTCAR','w')
    #orbcarfile=open(curdir+'/ORBCAR','w')
    wavecarfile=open(curdir+'/WAVECAR','w')


except:
    print("INCAR file not found\n")
    quit()

#default values
e=0.2
p=8
spin=0.8
M=1e6
tottime=2e6
d1=0
d2=0
d3=0

#read incar file
outcarfile.write('------------------------------------------\n')
outcarfile.write('Reading INCAR file...\n')
control=readINCAR(curdir+'/INCAR')
outcarfile.write('INCAR file read.\n\nConfiguration shown below:\n')

e=float(control.ecc)
p=float(control.p)
iota=float(control.iota)
spin=float(control.spin)
M=float(control.mass)
tottime=float(control.total_time)
THETA=float(control.theta)
PHI=float(control.phi)
mass_ratio=float(control.mass_ratio)
R_pc=float(control.r_pc)

if control.metric.lower() != 'kerr':
    defpar_val=np.array(control.defpar_val)
    if control.metric.lower() == 'krz' or control.metric.lower() == 'custom':
        d1=defpar_val[0]
        d2=defpar_val[1]
        d3=defpar_val[2]

outcarfile.write("ECC".ljust(16) +'%.10e'%e + '\n')
outcarfile.write("P".ljust(16) +'%.10e'%p + '\n')
outcarfile.write("IOTA".ljust(16) +'%.10e'%iota + '\n')
outcarfile.write("SPIN".ljust(16) +'%.10e'%spin + '\n')
outcarfile.write("MASS".ljust(16) +'%.10e'%M + '\n')
outcarfile.write("THETA".ljust(16) +'%.10e'%THETA + '\n')
outcarfile.write("PHI".ljust(16) +'%.10e'%PHI + '\n')
outcarfile.write("MASS_RATIO".ljust(16) +'%.10e'%mass_ratio + '\n')
outcarfile.write("R_PC".ljust(16) +'%.10e'%R_pc + '\n')
outcarfile.write("METRIC".ljust(16) + control.metric.upper() + '\n')

if control.metric.lower() != 'kerr':
    valstr='[ '
    for i in range(len(defpar_val)):
        valstr=valstr+'%.10e,'%(defpar_val[i])
    valstr=valstr+' ]'
    outcarfile.write("DEFPAR_VAL".ljust(16)+ valstr + '\n')

#orbit calculation
E,Lz,Q=getELQ(e,p,iota,spin)
outcarfile.write('\n\n')
outcarfile.write('------------------------------------------\n')
outcarfile.write('Orbit calculation start...\n')
#outcarfile.write('------------------------------------------\n')
print("METRIC: "+control.metric.lower())


os.system('cp '+XSPEGdir+'/trace '+curdir+'/a.out')
os.system('chmod 777 '+ curdir+'/a.out')    

os.system('echo ./a.out %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3,mass_ratio,e,p,iota,R_pc,THETA,PHI))

os.system('./a.out %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3,mass_ratio,e,p,iota,R_pc,THETA,PHI))
os.system('rm a.out')
if control.metric.lower() == 'custom':
    os.system('rm main.cpp')
    os.system('rm christoffel.cpp')
    os.system('rm def.h')
    os.system('rm def.h.gch')    

print('Orbit saved in '+curdir+'/ORBCAR')
#outcarfile.write('------------------------------------------\n')
outcarfile.write('Orbit Calculation done. Orbits saved in ORBCAR\n')
outcarfile.write('------------------------------------------\n')

#print('Computing waveform')
#waveform calculation
#outcarfile.write('\n\n')
#outcarfile.write('------------------------------------------\n')
#outcarfile.write('Waveform calculation start...\n')
#myt_sec,mywave=getwave('ORBCAR',THETA=THETA,PHI=PHI,M=M,mu=mass_ratio,R_pc=R_pc) #format: time(s) , h_plus + i * h_cross
#for ind in range(len(myt_sec)):
#    wavecarfile.write('%.10e \t%.10e \t%.10e \n'%(myt_sec[ind],np.real(mywave[ind]),np.imag(mywave[ind])))

outcarfile.write('Waveform calculation done. Waveform saved in WAVECAR\n')
outcarfile.write('------------------------------------------\n')
#print('Done')
print('Waveform saved in '+curdir+'/WAVECAR') 
