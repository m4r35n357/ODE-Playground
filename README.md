## News

### (26th October 2021)
To distinguish it from the MPFR code in master, the compendium branch has been renamed to pure_c: https://github.com/m4r35n357/ODE-Playground/tree/pure_c.
It still contains the Python (pure!) and Symplectic code, it is just a rename.

### (16th October 2021)
Python simulators and AD/TSM module have been removed to eliminate duplication.
They can now be found only on the compendium branch here: https://github.com/m4r35n357/ODE-Playground/tree/compendium, together with the long double version of TSM, Symplectic Integrator, and Black Hole code.

### (4th October 2021)

I have decided to update the MPFR implementation in line with the pure C version that I have been working on since I wrote the message below.
That remains in the compendium branch,and its README contains material not covered here.
The code now looks like a framework; model clients no longer need to provide the TSM algorithm, or any other conditional code.
The core TSM algorithm and the Taylor recurrences have all been brought up to date, and are now as minimal and readable as they can possibly be.

### Important (OBSOLETE - kept for historical reasons)
This branch is pretty much "finished", or "abandoned", current work (double precision, bifurcation diagrams, no function analysis) is happening  on the "compendium" branch here: https://github.com/m4r35n357/ODE-Playground/tree/compendium

### Background

This project is mainly a collection of programs in c and Python for evolving systems of ODEs using the Taylor Series Method (TSM), a rather old but poorly acknowledged technique based on forward mode Automatic Differentiation (AD).
TSM is a procedure for integrating ODEs using Taylor Series of arbitrary order, calculated to arbitrary precision (the former requires the latter in practice), using recurrence relations between time derivatives of increasing order.
It is therefore (or should be, if it were better known!) a serious competitor to the fourth-order RK4 for the vast majority of cases.

The code itself is tiny (as are the dynamic executables) and has been partly developed on, and is suitable for running on, a Raspberry Pi 4 computer.
The c code supports arbitrary precision (using the GNU MPFR library: https://www.mpfr.org), and is most suited to solving systems of ODEs to high accuracy.
The Python code uses float precision, and is most suited to interactive analysis and plotting of functions and their derivatives.

For the uninitiated, here is a review of the TSM itself and its history of repeated "discovery" and re-branding.
https://arxiv.org/abs/1111.7149

My work was inspired by parts of the following paper and its associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

Their software (called "Taylor") uses code generation to produce the recurrences and code to drive them, and it inspired me to try coding something more direct by hand.
These programs are the result.

Finally, more general resources on automatic differentiation can be found at the following portal: http://www.autodiff.org/ (see specifically the "Applications" and "Tools" sections).

## Taylor Series ODE Solvers, c/MPFR in master branch, pure c and Python in pure_c branch)

My primary aim was to be able to solve coupled nonlinear equations and investigate chaotic systems, without relying on "black-box" ODE solvers.
The resulting c code takes the form of a small (<200 loc) arbitrary precision Taylor Series "library", and the model-specific ODE simulators are tiny client programs to this, typically 25-35 loc each.
The header file taylor-ode.h contains a terse but complete description of the Taylor Series Method as implemented here, together with derivations of the recurrences that enable analysis of complex composed functions.

I have also duplicated the ODE solving functionality in Python 3 (at float precision), see the pure_c branch, but with extra testing and more advanced function analysis features (enabled by operator overloading and the Python REPL).

The recurrence rules (the "t-functions" in c and Python) are the key to calculating high order derivatives accurately, without needing finite differences.
They generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order values.
The functions provided cover the basic algebraic operations on Taylor Series (+ - * /), and also include several common functions:
* abs
* mul
* sqr
* sqrt
* div, inv
* exp
* sin(h)_cos(h)
* tan(h)_sec^2(h)
* pwr (f(x)^a, where a is real)
* ipwr (f(x)^a, where a is an integer)
* ln
* asin(h), acos(h), atan(h)

The recurrence relations used here are derived along the lines of (amongst other sources) http://www2.math.uni-wuppertal.de/wrswt/preprints/prep_05_4.pdf and http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).

There are also convenient factories for generating derivative "jets" of arbitrary order, and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, when properly called,  are all that is needed to solve systems of ODEs.
There is a fairly extensive collection of nonlinear ODE examples already implemented, in the file tsm-lorenz-dbg, and in the tsm-\*-*.c files.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and many others.

## Accuracy of solutions
Rather than estimating local error using the usual Taylor Series method, I have provided a clean numerical simulation (CNS) shell script, that makes it easy to run a "better" simulation alongside.
The differences can be plotted for easy visual comparison.

## Build/Test Environment (Debian/Ubuntu/Raspbian)

OS level requirements:
```
sudo apt install bc git build-essential pkg-config mesa-utils-extra python3-tk python3-dev libmpc-dev libfreetype6-dev libatlas-base-dev virtualenvwrapper gnuplot-x11 lcov
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground
cd ODE-Playground
```

#### Python 3 Packages (for plotting), please use a virtual environment!
There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python executables, which share output "formats").
In the example invocations given below, communication between the executable and plotting script uses a Unix pipe.
The dependencies are:
* matplotlib for 2D graphs and progressive graphs
* pi3d for 3D progressive trajectories
* gnuplot for static 3D trajectories

This is a typical Python virtual environment setup:
```
mkvirtualenv --python /usr/bin/python3 ad
pip install matplotlib pillow pi3d pytest pytest-cov mutmut ipython
```

#### Running Python Tests
Testing Requirements
* pytest
* pytest_cov
* mutmut (optional)

Most of the code is covered several times over.
Tests have been mutation tested with mutmut.
```
pytest ad_test.py solver_test.py -v
```

#### c Build (GCC or Clang)
```
./clean
./build [clang]
```
There should be NO errors or warnings.  [UPDATE: kerr-image.c shows warnings on arm64; it is 3rd party code]

#### Running c Tests
```
./libad-test-dbg 32 20 1 1e-18

Total: 33, PASSED 33
```
The final parameter can be set to 0 (or left absent) for a summary, 1 for individual tests, or 2 for full detail of Taylor Series.
Depending on the x value, some tests might be skipped owing to domain restrictions on some of the functions involved.

libad-test-dbg c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | (approximate) precision in bits
2 | order of Taylor Series
3 | x value
4 | tolerance
5 | debug level (optional)

#### Code coverage
```
./coverage
```
The output contains file system links to the HTML results

#### C Code profiling
```
./profile
```
The results are printed to stdout

#### Find examples for ODE parameters and other things:
```
grep Example *
```
## Solving and Plotting ODEs
This use case only involves calling the "t-functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm-lorenz-dbg and tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.
To find some example invocations:
```
grep Example *
```
Where CPU timings are given, they are made on a Raspberry Pi 400, mildly overclocked to 2100MHz, and writing output to a tmpfs file.

#### Run a basic ODE simulation (ODE call):

Runs a named simulation, and prints results to stdout.
Each output line consists of a column each for x, y, z, t, followed by cumulative CPU usage

**tsm-model-type** (c executables)

Parameter | Meaning
----------|-----------
1 | x,y,z display precision in bits
2 | (approximate) internal precision in decimal places
3 | order of Taylor Series
4 | time step
5 | number of steps
6,7,8 | initial conditions, x0,y0,z0
9+ | Model parameters

##### Run & plot (3D plot using pi3d):
```
./tsm-thomas-dbg 9 113 10 0.1 30000 1 0 0 .185 | ./plot3d.py
```
##### Run & plot (animated matplotlib graph):
```
./tsm-lorenz-dbg 9 113 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-dbg 9 113 10 0.1 30000 1 0 0 .185 >/tmp/$USER/data
gnuplot -p -e "set xyplane 0; set view 54.73561,135; set xlabel 'X'; set ylabel 'Y'; set zlabel 'Z'; splot '/tmp/ian/data' with lines"
```
##### Run & plot (2D gnuplot graph):
```
./tsm-lorenz-dbg 9 113 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/$USER/data
gnuplot -p -e "set terminal wxt size 1200,900; plot '/tmp/$USER/data' using 4:1 with lines, '/tmp/$USER/data' using 4:2 with lines, '/tmp/$USER/data' using 4:3 with lines"
```
It should be possible to send output directly to gnuplot via a pipe, but many versions segfault when reading stdin so I now specify a temporary file instead.

#### Clean Numerical Simulation:

In a chaotic system, accuracy can only be maintained for a finite simulation time.
This script runs a given simulation twice, the second time with a "better" integrator, and shows the differences graphically.
The maximum achievable clean simulation time is determined by machine precision alone.
Step size and order can affect the efficiency of the calculations, but not the ultimate accuracy!

**cns** (shell script)

Parameter | Meaning
----------|-----------
1 | CNS function, Selects a better integrator for comparison, see below
2 | deviation threshold
3+ | ODE call (you can now use precision "names" in cns scripts)

CNS function | Meaning
----------|-----------
step2 | The step size is halved (this is  now the _only_ "better" integrator!)

##### CNS plot (matplotlib diff graph):

```
./cns step2 1 ./tsm-lorenz-static 9 long 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
```
#### Example output - ~300 time units
```
$ time -p ./cns step2 1 ./tsm-lorenz-static 6 408 102 .01 32000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 408 102 .01 32000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 408 102 .005000 64000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 408 bits
 MPFR default precision: 408 bits
threshold: 1.0e+00  t: 301.330  cpu: 103.982 208.085
real 224.00
user 330.71
sys 0.81
```
#### Example output - ~450 time units
```
$ time -p ./cns step2 1 ./tsm-lorenz-static 6 630 153 .01 46000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 630 153 .01 46000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 630 153 .005000 92000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 630 bits
 MPFR default precision: 630 bits
threshold: 1.0e+00  t: 453.560  cpu: 545.447 1091.855
real 1109.91
user 1660.16
sys 0.72
```
#### Example output - ~600 time units
```
$ time -p ./cns step2 1 ./tsm-lorenz-static 6 840 204 .01 62000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 840 204 .01 62000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 840 204 .005000 124000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 840 bits
 MPFR default precision: 840 bits
threshold: 1.0e+00  t: 600.870  cpu: 1956.807 3917.292
real 4052.25
user 6059.11
sys 2.41
```

(matplotlib plot not shown!)

#### CNS Duration Scanning

Runs a simulation repeatedly, incrementing the order of integration, and for each order shows the simulation time when the deviation threshold is exceeded.
You can run this to determine the maximum _useful_ integrator order to use, for a given precision and step size.

**cns-scan** (shell script) 

Parameter | Meaning
----------|-----------
1 | minimum order for Taylor integrator
2 | maximum order for Taylor integrator
3 | deviation threshold
4+ | ODE call (the order parameter should be a placeholder, e.g. "_", to avoid confusion)

#### CNS duration vs. Simulation Order (gnuplot graph) for the given step size:

The following commands perform a scan, and plot the simulation time and cpu time as histograms against integrator order.
Quadruple precision can be clean up to ~76 time units:
```
./cns-scan 4 32 1 ./tsm-lorenz-static 6 quad _ .01 10000 -15.8 -17.48 35.64 10 28 8 3 | tee /tmp/$USER/data
```

Octuple precision can be clean up to ~173 time units:
```
./cns-scan 24 60 1 ./tsm-lorenz-static 6 oct _ .01 18000 -15.8 -17.48 35.64 10 28 8 3 | tee /tmp/$USER/data
```

To plot clean simulation time and CPU vs. order:
```
gnuplot -p -e "set ytics nomirror; set y2tics; plot '/tmp/$USER/data' using 1:2 axes x1y1 with boxes, '/tmp/$USER/data' using 1:3 axes x1y2 with boxes"
```

To plot required order and CPU time against (desired) maximum clean simulation time from the same data:
```
gnuplot -p -e "set ytics nomirror; set y2tics; plot '/tmp/$USER/data' using 2:1 axes x1y1 with points, '/tmp/$USER/data' using 2:3 axes x1y2 with points"
```

## ic script - Sensitivity to variation in initial conditions
This script is used to generate deviation data for chaos scanning, but the data can also be plotted in real time using matplotlib.
As well as the trajectory specified in the command arguments, six others are created and evolved; each one is the centre of the face of a cube around the original value

ic shell script ||

Parameter | Meaning
----------|-----------
1 | Initial separation between "original" trajectory and the extra ones
2+ | ODE call

The simulation is run seven times in parallel processes, the original along with each perturbed x, y, z.
```
./ic .001 ./tsm-lorenz-dbg 9 53 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 2>/dev/null
oo ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
x+ ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.799 -17.48 35.64 10 28 8 3
x- ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.801 -17.48 35.64 10 28 8 3
y+ ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.8 -17.479 35.64 10 28 8 3
y- ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.8 -17.481 35.64 10 28 8 3
z+ ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.8 -17.48 35.641 10 28 8 3
z- ./tsm-lorenz-dbg 9 9 53 10 .01 10001 -15.8 -17.48 35.639 10 28 8 3
```
(3D plot not shown)
