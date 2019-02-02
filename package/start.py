import numpy
import numpy as np
from KRZmetric import *
import os

XSPEGdir=os.environ['XSPEGLIB']

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
    if control.metric.lower == 'krz':
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
    outcarfile.write("DEFPAR_VAL".ljust(16)+ defpar_val + '\n')

#orbit calculation
E,Lz,Q=getELQ(e,p,iota,spin)
outcarfile.write('\n\n')
outcarfile.write('------------------------------------------\n')
outcarfile.write('Orbit calculation start...\n')
#outcarfile.write('------------------------------------------\n')

os.system('cp '+XSPEGdir+'/trace '+curdir  )
os.system('echo ./trace %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e \n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3))

os.system('./trace %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e \n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3))
os.system('rm trace')

print('Orbit saved in '+curdir+'/ORBCAR')
#outcarfile.write('------------------------------------------\n')
outcarfile.write('Orbit Calculation done. Orbits saved in ORBCAR\n')
outcarfile.write('------------------------------------------\n')

print('Computing waveform')
#waveform calculation
outcarfile.write('\n\n')
outcarfile.write('------------------------------------------\n')
outcarfile.write('Waveform calculation start...\n')
myt_sec,mywave=getwave('ORBCAR',THETA=THETA,PHI=PHI,M=M,mu=mass_ratio,R_pc=R_pc) #format: time(s) , h_plus + i * h_cross
for ind in range(len(myt_sec)):
    wavecarfile.write('%.10e \t%.10e \t%.10e \n'%(myt_sec[ind],np.real(mywave[ind]),np.imag(mywave[ind])))

outcarfile.write('Waveform calculation done. Waveform saved in WAVECAR\n')
outcarfile.write('------------------------------------------\n')
print('Done')
print('Waveform saved in '+curdir+'/WAVECAR') 
