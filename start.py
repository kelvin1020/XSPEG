# coding: utf-8
import numpy
from numpy import sin,cos,sqrt
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import scipy

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


# XSPEGdir=os.environ['XSPEGLIB']
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
