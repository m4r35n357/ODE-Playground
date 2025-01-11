## News

### (December 2025)

* (pure_c) Arbitrary-precision solver and models now using GNU bc (rather than previous MPFR in master branch).
* (pure_c) Animated ODE plots added (using Python/Matplotlib)
* (python) Python branch now just five models (using TSM & RK4), plus interactive shell and function plotters/solvers
* (master) MPFR branch much `updated
* Updated 'ode-playground' scripts

### (June 2022)
Since last update:

branch | contents
----------|-----------
master | MPFR/c99 code for TSM
pure_c | Pure c99 code for TSM & RK4, plus symplectic integration and dual numbers for Hamiltonians
pure_c_python | all the above, plus Python TSM/RK4, plotters, and interactive solvers/plotters for IPython

The arc/ar-functions have been added (and fully tested).

Documentation has been slightly improved.

Several cut & paste examples are provided in comments at the top of the following source files:
* tsm-lorenz.c
* tsm-thomas.c
* tsm-halvorsen.c

A blank initial character is provided to keep these long lines off your history (bash).
These examples allow simple selection of static & debug executables, and Python versions (where valid).

There are also yad dialogues for running tests and making executables using zig-build:
* libad-test.c


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

This project is mainly a collection of programs in c and Python (in the pure_c_python branch only!) for evolving systems of ODEs using the Taylor Series Method (TSM), a rather old but poorly acknowledged technique based on forward mode Automatic Differentiation (AD).
TSM is a procedure for integrating ODEs using Taylor Series of arbitrary order, calculated to arbitrary precision (the former requires the latter in practice), using recurrence relations between time derivatives of increasing order.
It is therefore (or should be, if it were better known!) a serious competitor to the fourth-order RK4 for the vast majority of cases.

The code itself is tiny (as are the dynamic executables) and has been partly developed on, and is suitable for running on, a Raspberry Pi 4 computer.
The c code supports arbitrary precision (using the GNU MPFR library: https://www.mpfr.org), and is most suited to solving systems of ODEs to high accuracy.
The Python code (now in the pure_c_python branch only) uses float precision, and is most suited to interactive analysis and plotting of functions and their derivatives.

For the uninitiated, here is a review of the TSM itself and its history of repeated "discovery" and re-branding.
https://arxiv.org/abs/1111.7149

My work was inspired by parts of the following paper and its associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

Their software (called "Taylor") uses code generation to produce the recurrences and code to drive them, and it inspired me to try coding something more direct by hand.
These programs are the result.

Finally, more general resources on automatic differentiation can be found at the following portal: http://www.autodiff.org/ (see specifically the "Applications" and "Tools" sections).

## Taylor Series ODE Solvers, c/MPFR in master branch, pure c and Python in pure_c_python branch)

My primary aim was to be able to solve coupled nonlinear equations and investigate chaotic systems, without relying on "black-box" ODE solvers.
The resulting c code takes the form of a small (<200 loc) arbitrary precision Taylor Series "library", and the model-specific ODE simulators are tiny client programs to this, typically 25-35 loc each.
The header file taylor-ode.h contains a terse but complete description of the Taylor Series Method as implemented here, together with derivations of the recurrences that enable analysis of complex composed functions.

I have also duplicated the ODE solving functionality in Python 3 (at float precision), see the pure_c_python branch, but with extra testing and more advanced function analysis features (enabled by operator overloading and the Python REPL).

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

## Quick Start

### Requirements - Debian/Ubuntu packages
```
sudo apt install bc git build-essential libmpc-dev libfreetype6-dev gnuplot-qt lcov
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground
cd ODE-Playground
```

#### c Build (GCC by default, Clang optional)
```
./clean
./build [clang]
```
There should be NO errors or warnings.  [UPDATE: kerr-image.c shows warnings on arm64; it is 3rd party code]

[UPDATE] clang shows warnings for the unsupported -Wunsuffixed-float-constants warning!

Each client is built _both_ as a dynamic executable with asserts and debug symbols, and as a stripped static executable with asserts disabled.

#### Running c Tests

The tests enforce a redundant web of densely interrelated functionality that cannot exist in the presence of coding errors! ;)

**libad-test** (c executable)

Parameter | Meaning
----------|-----------
1 | Display precision
2 | Precision in bits
3 | Order
4 | X value
5 | Error limit (absolute)
6 | (Optional) verbosity: 0: summary (default), 1: list, 2: detail

The final parameter can be set to 0 (or left absent) for a summary, 1 for individual tests, or 2 for full detail of Taylor Series.
Depending on the x value, some tests might be skipped owing to domain restrictions on some of the functions involved.

```
$ ./libad-test 15 237 64 -.5 1e-64
argc: 6, argv: [ ./libad-test 15 237 64 -.5 1e-64 ]
Taylor Series Method 
+1.000000000000000e+00 +1.000000000000000e+00 +1.000000000000000e+00 +0.000000000e+00 0.000
+1.105170918075648e+00 +1.000000000000000e+00 +9.048374180359596e-01 +1.000000000e-01 0.000
+1.221402758160170e+00 +1.000000000000000e+00 +8.187307530779819e-01 +2.000000000e-01 0.000
+1.349858807576003e+00 +1.000000000000000e+00 +7.408182206817179e-01 +3.000000000e-01 0.000
+1.491824697641270e+00 +1.000000000000000e+00 +6.703200460356393e-01 +4.000000000e-01 0.000
+1.648721270700128e+00 +1.000000000000000e+00 +6.065306597126334e-01 +5.000000000e-01 0.001
+1.822118800390509e+00 +1.000000000000000e+00 +5.488116360940264e-01 +6.000000000e-01 0.001
+2.013752707470477e+00 +1.000000000000000e+00 +4.965853037914095e-01 +7.000000000e-01 0.001
+2.225540928492468e+00 +1.000000000000000e+00 +4.493289641172216e-01 +8.000000000e-01 0.001
+2.459603111156950e+00 +1.000000000000000e+00 +4.065696597405991e-01 +9.000000000e-01 0.001
+2.718281828459045e+00 +1.000000000000000e+00 +3.678794411714423e-01 +1.000000000e+00 0.001
Check: e^1  e^0  e^-1
+2.718281828459045e+00 +1.000000000000000e+00 +3.678794411714423e-01 +1.000000000e+00 0.002
Horner Summation 
 23   23.000
153  153.000
201  201.000
Recurrence Relations: x = -0.5
Total: 49, PASSED 28, SKIPPED 21
```

#### Code coverage
Creates a web page summary
```
make clean && make CCC=cov && make coverage
```
The output contains file system links to the HTML results

#### C Code profiling
Very basic information, included just for completeness
```
./profile gpt ./tsm-lorenz-std  6 237 16 .01 1000000  -15.8 -17.48 35.64  10 28 8 3
./profile gcc ./tsm-lorenz-std  6 237 16 .01 1000000  -15.8 -17.48 35.64  10 28 8 3
```
The results are printed to stdout

#### Find examples for ODE parameters and many other things:
Useful commands are frequently added to the comments in source headings.
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
Each output line consists of a column each for x, y, z, t, followed by three empty turning point tags and cumulative CPU usage

**tsm-model-type** (c executables)

Parameter | Meaning
----------|-----------
1 | x,y,z display precision in decimal places
2 | internal precision in bits
3 | order of Taylor Series
4 | time step
5 | number of steps
6,7,8 | initial conditions, x0,y0,z0
9+ | Model parameters

##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-dbg 6 113 10 0.1 30000 1 0 0 .185 >/tmp/$USER/data

 gnuplot << EOF
set xyplane 0
set view 54.73561,135
set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'
splot '/tmp/$USER/data' with lines
pause mouse close
print "Done"
EOF
```
##### Run & plot (2D gnuplot graph):
```
./tsm-lorenz-dbg 6 113 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/$USER/data

 gnuplot << EOF
set terminal wxt size 1200,900
plot '/tmp/$USER/data' using 4:1 with lines, '/tmp/$USER/data' using 4:2 with lines, '/tmp/$USER/data' using 4:3 with lines
pause mouse close
print "Done"
EOF
```
It should be possible to send output directly to gnuplot via a pipe, but many versions segfault when reading stdin so I now specify a temporary file instead.

### Bifurcation Diagrams:

See pure_c branch.

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
3+ | ODE call (you can now use precision "names" in CNS scripts)

CNS function | Meaning
----------|-----------
step2 | The step size is halved (this is  now the _only_ "better" integrator!)
nosim | User-defined comparison between /tmp/$USER/dataA and /tmp/$USER/dataB

NOTE: the CPU time required for a clean simulation is of order:
```
O(required precision * order^2 * clean simulation time)
```
Since (empirically, for these Lorenz parameters and timestep)
* clean simulation time is approximately proportional to the required precision
* clean simulation time is approximately proportional to the Taylor series order

the CPU requirement can also be seen as:
```
O(required precision^4) or O(order^4) or O(clean simulation time^4)
```
##### Make a CNS plot:

#### Example output - quadruple precision
```
./cns step2 1 ./tsm-lorenz-static 6 quad 28 .01 10000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 quad 28 .01 10000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 113 28 .005000 20000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 113 bits
 MPFR default precision: 113 bits
threshold: 1.0e+00  t: 76.440  cpu: 1.264 2.590
real 3.41
user 4.97
sys 0.10
```
#### Example output - octuple precision
```
time -p ./cns step2 1 ./tsm-lorenz-static 6 oct 58 .01 18000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 oct 58 .01 18000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 237 58 .005000 36000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 237 bits
 MPFR default precision: 237 bits
threshold: 1.0e+00  t: 171.910  cpu: 12.511 24.808
real 26.11
user 39.03
sys 0.17
```
#### Example output - ~300 time units
```
time -p ./cns step2 1 ./tsm-lorenz-static 6 408 102 .01 32000 -15.8 -17.48 35.64 10 28 8 3
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
time -p ./cns step2 1 ./tsm-lorenz-static 6 630 153 .01 46000 -15.8 -17.48 35.64 10 28 8 3
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
time -p ./cns step2 1 ./tsm-lorenz-static 6 840 204 .01 62000 -15.8 -17.48 35.64 10 28 8 3
Clean Numerical Simulation: [step2 1 ./tsm-lorenz-static 6 840 204 .01 62000 -15.8 -17.48 35.64 10 28 8 3]
Better: ./tsm-lorenz-static 6 840 204 .005000 124000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 840 bits
 MPFR default precision: 840 bits
threshold: 1.0e+00  t: 600.870  cpu: 1956.807 3917.292
real 4052.25
user 6059.11
sys 2.41
```
If you need to re-plot after closing gnuplot, either use the "nosim" argument, or:
```
 gnuplot << EOF                                                             
set key horizontal left
plot '/tmp/$USER/dataA' using 4:1 t 'xA' with lines lc black, '' u 4:2 t 'yA' w l lc black, '' u 4:3 t 'zA' w l lc black, '/tmp/$USER/dataB' using 4:1 t 'xB' with lines, '' u 4:2 t 'yB' w l, '' u 4:3 t 'zB' w l
pause mouse close
print "Done"
EOF
```

#### CNS Duration Scanning

Runs a simulation repeatedly, incrementing the order of integration, and for each order shows the simulation time when the deviation threshold is exceeded.
You can run this to determine the maximum _useful_ integrator order to use, for a given precision and step size.

**cns-scan** (shell script) 

Parameter | Meaning
----------|-----------
1 | Minimum order for Taylor integrator
2 | Maximum order for Taylor integrator
3 | deviation threshold
4+ | ODE call (with "order" argument set to "_", to avoid confusion)

#### CNS duration vs. Simulation Order (gnuplot graph) for the given step size:

The following commands perform a scan, and plot the simulation time and cpu time as histograms against integrator order.
Quadruple precision can be clean up to ~76 time units for Lorenz:
```
./cns-scan 4 32 1 ./tsm-lorenz-static 6 quad _ .01 10000 -15.8 -17.48 35.64 10 28 8 3 | tee /tmp/$USER/data
```

Octuple precision can be clean up to ~173 time units for Lorenz:
```
./cns-scan 24 60 1 ./tsm-lorenz-static 6 oct _ .01 18000 -15.8 -17.48 35.64 10 28 8 3 | tee /tmp/$USER/data
```

To plot clean simulation time and CPU vs. order:
```
 gnuplot << EOF
set key left
set ytics nomirror
set y2tics
set xlabel 'Taylor Series Order'
set ylabel 'CNS Time, model units'
set y2label 'CPU Time, seconds'
plot '/tmp/$USER/data' using 1:2 axes x1y1 title 'CNS' with boxes, '/tmp/$USER/data' using 1:3 axes x1y2 title 'CPU' with boxes
pause mouse close
print "Done"
EOF
```
