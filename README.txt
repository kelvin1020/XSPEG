XSPEG v0.1.0 

This software (as well as its future releases) will be  used to compute gravitational waveforms from Extreme-Mass-Ratio Inspirals (EMRIs) in non-Kerr spacetime, although currently the first version only inplement one kind of non-Kerr parametrization and radiation reaction is not inplemented (so no inspirals actually happen in v0.1.0, the orbits are geodesics).
The usage of this software will be very similar to VASP, a software widely used in solid state physics computation.

INSTALL:
run the routine "install_xspeg.sh" to install this software:
./install_xspeg.sh

USE:
1. prepare the INCAR file where you specify the parameters in the format explained in manual.pdf (see sample_input/ for an example).
2. run the software by the command:
xspeg
3. Three files will be generated in your current directory:
OUTCAR: explaining the parameters and the computing process
ORBCAR: the orbit of the test particle, in the format explained in manual.pdf (see sample_input/ORBCAR for an example).
WAVECAR: the gravitational waveform, in the format explained in manual.pdf (see sample_input/WAVECAR for an example).


Author: Xin Shuo, shuoxintj@gmail.com
