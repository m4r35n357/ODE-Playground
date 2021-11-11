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

## Taylor Series ODE Solvers (c/MPFR and Python)

My primary aim was to be able to solve coupled nonlinear equations and investigate chaotic systems, without relying on "black-box" ODE solvers.
The resulting c code takes the form of a small (<200 loc) arbitrary precision Taylor Series "library", and the model-specific ODE simulators are tiny client programs to this, typically 25-35 loc each.
The header file taylor-ode.h contains a terse but complete description of the Taylor Series Method as implemented here, together with derivations of the recurrences that enable analysis of complex composed functions.

I have also duplicated the ODE solving functionality in Python 3 (at float precision), but with extra testing and more advanced function analysis features (enabled by operator overloading and the Python REPL).

The recurrence rules (the "t-functions" in c and Python) are the key to calculating high order derivatives accurately, without needing finite differences.
They generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order values.
The functions provided cover the basic algebraic operations on Taylor Series (+ - * /), and also include several common functions:
* abs
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, where a is a scalar)
* ln
* asin(h), acos(h), atan(h) - Python only

The recurrence relations used here are derived along the lines of (amongst other sources) http://www2.math.uni-wuppertal.de/wrswt/preprints/prep_05_4.pdf and http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).

There are also convenient factories for generating derivative "jets" of arbitrary order, and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, when properly called,  are all that is needed to solve systems of ODEs.
There is a fairly extensive collection of nonlinear ODE examples already implemented, in the file tsm-lorenz-dbg, and in the tsm-\*-*.c files.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and many others.

## Accuracy of solutions
Rather than estimating local error using the usual Taylor Series method, I have provided a clean numerical simulation (CNS) shell script, that makes it easy to run a "better" simulation alongside.
The differences can be plotted for easy visual comparison.

## Function Analysis (Python)

Partly to verify my own implementation of the Taylor recurrence rules, I have added a demonstration of using series arithmetic to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Of course it is more generally useful for finding inverse values (where real solutions exist) of complicated functions, not just their roots.

The "higher-level" ad_functions (or the Series class in Python) manipulate entire  Taylor series "jets" at once, so are only useful for univariate functions.
In Python there is an overloaded operator "~" which extracts the actual derivative values from the taylor series coefficents, as well as additional operators for (negation and **).
The \*\* (power) operator caters for f(x)^a, a^f(x) and f1(x)^f2(x), subject to domain limitations on f(x).
There are also functions for (matching the t_functions):
* abs
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, where a is a scalar)
* ln
* asin(h), acos(h), atan(h) - Python only

Using these "higher level" functions, Newton's method is implemented trivially, but I have also provided an implementation of the bisection method for comparison.

The plotters.py script enables the analysis of a "model" function's derivatives along with the value itself, across a range of the input variable.
That file contains a selection of example invocations in the comments.
Optionally it will analyse that range for roots, extrema and inflection points using the lower order derivatives.
Here is a seventh-degree polynomial model by way of example, and a trigonometric identity as a check of the sin and sqr functions:
```python
def septic(a):
    return (a + 7) * (5 + a) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6)

def trig(a):
    return (3 * a).sin - 3 * a.sin + 4 * a.sin**3
```
Operators on or between jets and numbers are implemented by overloading, and the functions are implemented as _properties_ of a jet.
These two approaches produce a fairly readable format for coding the models.
The value parameter allows simple calculation of function inverse values, as well as roots.

In summary, there are five main areas of application for the code:
* solving nonlinear ODEs (to arbitrary precision with c/MPFR)
* scanning ODE parameters for detecting chaos in solutions
* checking global accuracy of solutions
* plotting functions and their (higher) derivatives, with solution, turning-point, and inflection analysis (Python is best here)
* various interactive investigations in the Python console

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
$ pytest ad_test.py solver_test.py -v
```

#### c Build (GCC or Clang)
```
$ ./clean
$ ./build [clang]
```
There should be NO errors or warnings.  [UPDATE: kerr-image.c shows warnings on arm64; it is 3rd party code]

#### Running c Tests
```
$ ./libad-test-dbg 32 20 1 1e-18

Total: 33, PASSED 33
```
The final parameter can be set to 0 (or left absent) for a summary, 1 for individual tests, or 2 for full detail of Taylor Series.
Depending on the x value, some tests might be skipped owing to domain restrictions on some of the functions involved.

libad-test-dbg c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | (approximate) precision in decimal places
2 | order of Taylor Series
3 | x value
4 | tolerance
5 | debug level (optional)

#### Code coverage
```
$ ./coverage
```
The output contains file system links to the HTML results

#### C Code profiling
```
$ ./profile
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
$ grep Example *.py tsm-*.c
```
#### Run an ODE simulation (ODE call):

Runs a named simulation, and prints results to stdout

tsm-\*-\* c executables ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | x, y, z output precision in decimal places (0 for full)
2 | (approximate) internal precision in decimal places
3 | order of Taylor Series
4 | time step size
5 | number of steps
6,7, 8 | initial conditions, x0, y0, z0
9+ | ODE parameters

##### Run & plot (3D plot using pi3d):
```
./tsm-thomas-dbg 9 32 10 0.1 30000 1 0 0 .185 | ./plot3d.py
```
##### Run & plot (animated matplotlib graph):
```
./tsm-lorenz-dbg 9 32 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
##### Run & plot (3D gnuplot graph):
```
./tsm-thomas-dbg 9 32 10 0.1 30000 1 0 0 .185 | gnuplot -p -e "set terminal wxt size 1200,900; splot '<cat' with lines"
```
##### Run & plot (2D gnuplot graph):
```
./tsm-lorenz-dbg 9 32 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/$USER/data
gnuplot -p -e "set terminal wxt size 1200,900; plot '/tmp/$USER/data' using 4:1 with lines, '/tmp/$USER/data' using 4:2 with lines, '/tmp/$USER/data' using 4:3 with lines"
```

#### Clean Numerical Simulation:

This is a relatively new approach to dealing with the global error of ODE simulations, described in detail here: https://arxiv.org/abs/1109.0130.
As a rough guide to the accuracy of a solution, it can be compared with a "better" solution (one made with "better" solver parameters), and discarded at the point where the solutions diverge.
The two simulations are run in parallel processes, but obviously the "better" solution takes longer.
There are three alternative strategies for creating the "better" integrator:

cns shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | CNS function, Selects a better integrator for comparison, see below
2+ | ODE call

CNS function Parameter | Meaning
----------|-----------
step2 | The step size is halved
both | The order is increased by one, and the step size by one half

##### CNS plot (matplotlib diff graph):
```
./cns both ./tsm-lorenz-static 9 32 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
```
#### Example output - 300 time units
```
$ ./cns step2 ./tsm-lorenz-dbg 15 130 102 .01 35000 -15.8 -17.48 35.64 10 28 8 3
Better: ./tsm-lorenz-dbg 138 103 .005000 70000 -15.8 -17.48 35.64 10 28 8 3
 MPFR default precision: 458 bits
 MPFR default precision: 431 bits
Threshold: 1.0e-12, t: 267.060
Threshold: 1.0e-09, t: 278.700
Threshold: 1.0e-06, t: 284.950
Threshold: 1.0e-03, t: 293.140
Threshold: 1.0e+00, t: 301.320
```
(matplotlib plot not shown!)

600 time units
```
$ ./cns step2 ./tsm-lorenz-dbg 15 240 204 .01 65000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

1500 time units
```
$ ./cns step2 ./tsm-lorenz-dbg 15 800 501 .005 150000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

#### CNS Duration Scanning

Runs a simulation repeatedly with increasing order of integration, for each order showing the simulation time when the deviation threshold is exceeded.

cns-scan (shell script) ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | Precision, also used as maximum order for Taylor integrator (minimum is 1)
2 | deviation threshold
3+ | ODE call (the precision and order parameters should be placeholders to avoid confusion)

##### CNS duration vs. Simulation Order (gnuplot graph):
```
$ ./cns-scan step2 32 1 ./tsm-lorenz-static 6 _ _ .01 10000 -15.8 -17.48 35.64 10 28 8 3  | gnuplot -p -e "plot '<cat' with boxes"
```

## ic script - Sensitivity to variation in initial conditions
This script is used to generate deviation data for chaos scanning, but the data can also be plotted in real time using matplotlib.
As well as the trajectory specified in the command arguments, six others are created and evolved; each one is the centre of the face of a cube around the original value

ic shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
separation | Initial separation between "original" trajectory and the additional ones

The simulation is run seven times in parallel processes, the original along with each perturbed x, y, z.
```
$ ./ic .001 ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
oo ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
x+ ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.799 -17.48 35.64 10 28 8 3
x- ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.801 -17.48 35.64 10 28 8 3
y+ ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.479 35.64 10 28 8 3
y- ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.481 35.64 10 28 8 3
z+ ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.48 35.641 10 28 8 3
z- ./tsm-lorenz-dbg 9 9 32 10 .01 10001 -15.8 -17.48 35.639 10 28 8 3
```
(3D plot not shown)
