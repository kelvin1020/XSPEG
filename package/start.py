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

if control.metric.lower() == 'krz'or control.metric.lower() == 'kerr':
    #print("METRIC: "+control.metric.upper)
    os.system('cp '+XSPEGdir+'/trace '+curdir  )
elif control.metric.lower() == 'custom':
    os.system('cp '+ XSPEGdir+'/custom/main.cpp ' + curdir)
    os.system('cp '+ XSPEGdir+'/custom/def.h ' + curdir)
    os.system('cp '+ XSPEGdir+'/custom/christoffel.cpp ' + curdir)
    #os.system('echo cp')
    #os.system('echo cp '+ XSPEGdir+'/custom/main.cpp ' + curdir)
    #os.system('echo cp '+ XSPEGdir+'/custom/def.h ' + curdir)
    #os.system('echo cp '+ XSPEGdir+'/custom/christoffel.cpp ' + curdir)

    try:
        os.system('g++ '+curdir+'/*.cpp '+curdir+'/*.h -o trace')
        #os.system('echo g++ '+curdir+'/*.cpp '+curdir+'/*.h -o trace')

    except:
        print('Please check custom metric files you prepared')
os.system('echo ./trace %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3,mass_ratio,e,p,iota))

os.system('./trace %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n'%(M,spin,E,Lz,Q,p/(1-e),tottime,d1,d2,d3,mass_ratio,e,p,iota))
os.system('rm trace')
if control.metric.lower() == 'custom':
    os.system('rm main.cpp')
    os.system('rm christoffel.cpp')
    os.system('rm def.h')

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
