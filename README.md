XSPEG v0.3.0-limited
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kelvin1020/XSPEG/master)

This software will be  used to compute gravitational waveforms from Extreme-Mass-Ratio Inspirals (EMRIs) in general stationary, time independent spacetime.
"KRZ" metric is inplemented in current version. For KRZ metric, see reference [1] and [2].
2PN approximation with some corrections is adopted when evaluating radiative force, see reference [3] and manual.pdf.
The usage of this software will be very similar to VASP, a software widely used in solid state physics computation.

RUN:
run the script "start.py" to start this software:

python start.py

USE:
1. prepare the INCAR file where you specify the parameters in the format explained in manual.pdf.
*Custom metric are not supported in this edtion (only built-in metrics are available,  i.e. KRZ metric and Kerr metric)

2. Three files will be generated in your current directory:
OUTCAR: explaining the parameters and the computing process
ORBCAR: the orbit of the test particle, in the format explained in manual.pdf (see sample_output/ORBCAR for an example).
WAVECAR: the gravitational waveform, in the format explained in manual.pdf (see sample_output/WAVECAR for an example).

[1] Roman Konoplya, Luciano Rezzolla, and Alexander Zhidenko. General parametrization of axisymmetric black holes in metric theories of gravity. Phys. Rev. D, 93:064015, Mar 2016.
[2] Yueying Ni, Jiachen Jiang, and Cosimo Bambi. Testing the kerr metric with the iron line and the krz parametrization. Journal of Cosmology and Astroparticle Physics, 2016(09):014, 2016.
[3] Stanislav Babak, Hua Fang, Jonathan R. Gair, Kostas Glampedakis, and Scott A. Hughes. ¡°kludge¡± gravitational waveforms for a test-body orbiting a kerr black hole. Phys. Rev. D, 75:024005, Jan 2007.

Author: Xin Shuo, shuoxintj@gmail.com
