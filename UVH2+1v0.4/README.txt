Causal Viscous Hydro Code for Non-Central Heavy Ion Collisions

Authors:

Ulrike Romatschke, Paul Romatschke and Matthew Luzum

version 0.2 February 2009
version 0.3 April 2010
version 0.3 September 2010


----------------------------------------------------------------

Compiling the code:


* the hydro-code is in the zipped-tar file named 

UVH2+1.tgz

It uses the gsl-library obtainable from 

http://www.gnu.org/software/gsl/

(versions newer than 1.8 are now OK wih version 0.3)

The hydro code expects the link "gsl" in the directory to point to 
the gsl directory.

You should then be able to compile the different modules of the code:

make initE
make vh2
make convert
make extract 

* the particle decay routines (adapted from AZHYDRO) are in the
zipped-tar file named 

reso.tgz

Put it into a separate directory where you should be able to compile
it with:

make reso 

Then take the executable "reso" and put it into the hydro directory .

* the CGC initial conditions (adapted from Y. Nara's fKLN code)
are in the zipped-tar file named

kln-2.21forUVH2+1v0.4.tgz

Put it into a separate subdirectory of the hydro code (e.g. in the
subdirectory "kln") and compile with "make". Put the executable
"kln.exe" into the hydro directory.


---------------------------------------------------------------

Using the code:

* Read the disclaimer in UVH2+1.cpp

* The code is separated into several modules:

 1a) initE: this module produces the Glauber initial energy-density
distribution in the transverse plane. The current version
lets you play with different versions of the Glauber initial
condition parametrization. 

Input: data/params.txt
Output: initE produces the file "inited.dat"

1b) kln.exe: this module (adapted from Y.Nara's fKLN code)
produces the CGC initial energy-density
distribution in the transverse plane. Note that this module produces a lot
of output, including some warnings/suggestions, which we believe
can be ignored in most cases.

Input: params.txt, qcdIEOS.dat
Main Output: kln.exe produces the file "inited.dat" for use with vh2
Extra Output: "initgp.dat" in a format for use with gnuplot splot

2) vh2: this is the main hydro code. It takes the output
from running initE as an input for the energy density
and then solves the hydro equations until freeze-out is finished.
It puts snapshots of things like the energy-density profile
into the directory "data/snapshot" (if present--takes up much more
disk space this way), and the eccentricities into "ecc.dat" in "data/"

Input: params.txt, inited.dat
Main Output: freezeout.dat
Extra Output: ecc.dat

3) convert (for isothermal freezeout): this module converts the hydrodynamic 
degrees of freedom (energy density, fluid velocities,...) 
into particle spectra using the Cooper-Frye freeze-out prescription

Input: params.txt, freezeout.dat, pasin.dat, pasim.dat, pasinames.dat, gslist.dat
Output: phipspectra.dat, ptarr.dat

3a) prereso: (optional) same as "extract" below, but before resonance feed down

Input: params.txt, phipspectra.dat
Output: data/results/prereso*

4) reso: this module (adapted from the AZHYDRO code) takes the
particle spectra and calculates decays of unstable particles
using the particle data book values.

Input: reso.inp, freezeout.dat
Output: spec_ and PT_ files in "data/results"
(see pdg_weak.dat for a clue what the numbers mean)

5) extract: this module finally takes the stable particle
spectra and calculates dN/dy/d^2pt (alias "v0"), "v2" and "v4"
from them. These can then be further manipulated to give
the total multiplicity and integrated/min-bias v2

Output: data/results/*sv?.dat

* All freezeout routines (convert, prereso, reso, extract) 
calculate spectra in only one quadrant of the transverse
plane, appropriate when there is reflection symmetry across
the x and y axes.  Verson 0.4 introduces new routines that 
integrate over the entire transverse plane, for when asymmetric
initial conditions are used, denoted by appending "full" to the 
routine name (convertfull, preresofull, resofull, and extractfull).
These routines take 4 times as long to run.

* The main control-file of the code is called "params.txt"
and can be found in the "data/" subdirectory. It has
most of the adjustable switches in there and a minimal description
to them. The physics parameters should be self-explanatory anyway.

* For convenience, once you've specified the relevant
parameters you're interested in in "params.txt", you can
simply invoke the shell-script "exshell", which should
do the steps 1)-5) for you and you can simply look at the
results in "data/results".

----------------------------------------------------------------

If something does not work or you think that we forgot to
mention something important, let us know!
Have fun.


