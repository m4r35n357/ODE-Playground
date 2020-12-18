## Important
This branch is pretty much "finished", or "abandoned", current work (double precision, bifurcation diagrams, no function analysis) is happening  on the tiny branch here: https://github.com/m4r35n357/ODE-Playground/tree/tiny

## Background

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
There is a fairly extensive collection of nonlinear ODE examples already implemented, in the file tsm.py, and in the tsm-\*-*.c files.
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

#### Graph Plotting
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
Now you can just use it "in place" in your virtual environment or, optionally, build and install an ad package to your venv
```
python3 setup.py sdist
cd dist
tar xzvf ad-1.0.tar.gz
cd ad-1.0
python3 setup.py install
cd ../..
```

#### Running Python Tests
Testing Requirements
* pytest
* pytest_cov
* mutmut (optional)

Most of the code is covered several times over.
Tests have been mutation tested with mutmut.
```
$ pytest --cov=ad --cov-report html:cov_html ad_test.py solver_test.py -q
.......................................................................... [ 24%]
.......................................................................... [ 48%]
.......................................................................... [ 73%]
.......................................................................... [ 97%]
........                                                                   [100%]

```

#### c Build (GCC or Clang)
```
$ ./build
$ ./build clang
```

#### Running c Tests
##### Newer tests
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

##### Older tests
```
$ ./ad-test-dbg 6 2 1
```
Output originally designed for visual checking, or see below how to do it automatically (hard-coded to quad precision).

ad-test-dbg c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | order of Taylor Series
2 | x value
3 | y value

##### Big build and test command (recommended):
```
time -p ./build && ./ad-test-dbg 6 2 1 >/tmp/ad-test.txt; diff --context=1 /tmp/ad-test.txt ad-test.txt && ./libad-test-dbg 32 20 1 1e-18 1 && echo OK
```
There is also a sample Git pre-commit script that you might want to copy (and edit) to .git/hooks

##### c Code Coverage
```
rm -f *.gcno *.gcda
./build
./libad-test-dbg 9 32 20 2 1 1e-18
lcov --capture --directory . --output-file coverage.info
genhtml coverage.info --output-directory out
```
Then run tests and/or programs to generate data files, HTML at out/index.html

## Solving and Plotting ODEs
This use case only involves calling the "t-functions" in ad.py or taylor-ode.c.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm.py and tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.
To find some example invocations:
```
$ grep Example *.py tsm-*.c
```
Matplotlib progressive graph plotting in Python (second parameter ignored as python floats are fixed double precision):
```
$ ./tsm.py lorenz 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
or in c (notice absence of the model parameter and presence of the internal precision arg):
```
$ ./tsm-lorenz-dbg 15 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```

#### tsm.py and tsm-\*-\* Parameter Reference
The Python ODE solver, tsm.py, comprises a long "if" statement containing a "zoo" of pre-programmed ODE systems.
tsm-\*-\* represents the individual c solvers.

##### Python (float precision)
tsm.py ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | model (ODE) name
2 | ignored
3 | order of Taylor Series
4 | step size
5 | number of steps
6,7,8 | x0, y0, z0
9+ | ODE parameters

##### c (MPFR arbitrary precision)
tsm-\*-\* c executable ||
----------|-----------
Parameter | Meaning
----------|-----------
1 | x, y, z output precision in decimal places (0 for full)
2 | (approximate) internal precision in decimal places
3 | order of Taylor Series (plot interval in RK4)
4 | step size
5 | number of steps
6,7, 8 | x0, y0, z0
9+ | ODE parameters

Since the RK4 Lorentz simulator (c only) is by definition fixed order, the "order" parameter is used to pass in a "plot interval" (e.g. 10 means plot only every 10th result).

#### 3D static trajectory plotting (gnuplot)

Write to a data file
```
$ ./tsm.py lorenz 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/data
```
then plot it
```
$ gnuplot -p -e "splot '/tmp/data' with lines"
```

#### 3D progressive trajectory plotting (pi3d)
```
$ ./tsm.py lorenz 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plot3d.py
```

## cns script - Clean Numerical Simulation

This is a relatively new approach to dealing with the global error of ODE simulations, described in detail here: https://arxiv.org/abs/1109.0130.
As a rough guide to the accuracy of a solution, it can be compared with a "better" solution (one made with "better" solver parameters), and discarded at the point where the solutions diverge.
The two simulations are run in parallel processes, but obviously the "better" solution takes longer.
There are three alternative strategies for creating the "better" integrator:

cns shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
step | The step size is reduced by a quarter (suitable for RK4) comparisons
order | The Taylor Series order is increased by two
both | The order is increased by one, and the step size by one half

#### Example output - 300 time units
```
$ ./cns both ./tsm-lorenz-dbg 15 130 102 .01 35000 -15.8 -17.48 35.64 10 28 8 3
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
$ ./cns both ./tsm-lorenz-dbg 15 240 204 .01 65000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

1500 time units
```
$ ./cns both ./tsm-lorenz-dbg 15 800 501 .005 150000 -15.8 -17.48 35.64 10 28 8 3
```
(output not shown)

## ic script - Sensitivity to variation in initial conditions
This script is used to generate deviation data for chaos scanning, but the data can also be plotted in real time using matplotlib.
As well as the trajectory specified in the command arguments, six others are created and evolved; each one is the centre of the face of a cube around the original value

ic shell script ||
----------|-----------
Parameter | Meaning
----------|-----------
separation | Initial separation between "original" trajectory and the additional ones
plot / noplot | Allows non-interactive use without plotting the data (e.g. by the chaos scanner).

The simulation is run seven times in parallel processes, the original along with each perturbed x, y, z.
```
$ ./ic .001 plot ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
oo ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
x+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.799 -17.48 35.64 10 28 8 3
x- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.801 -17.48 35.64 10 28 8 3
y+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.479 35.64 10 28 8 3
y- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.481 35.64 10 28 8 3
z+ ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.641 10 28 8 3
z- ./tsm-lorenz-dbg 9 16 10 .01 10001 -15.8 -17.48 35.639 10 28 8 3
```
(3D plot not shown)

## Interactive Function Analysis with Python
This use case involves calling the Series methods and operators from ad.py, and the functions from plotters.py, from a Python interpreter.
This is how I like to set things up:
```
$ echo '
from math import *
from ad import *
from plotters import *
' > ipython.py
$ ipython -i ipython.py
```

##### Higher level function analysis
Plotting function and derivatives, together with root and turning point analysis:
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: f = lambda a: (a.exp + (a.sqr - 4.0).exp).ln

In [2]: scan_s(f)
NT  x: -1.962e+00  δx: +4.123e-16  f: +1.332e-15  \ ROOT___ 4
NT  x: -1.312e+00  δx: -0.000e+00  f: +0.000e+00  / MIN_MAX 4
NT  x: -1.849e-02  δx: -2.714e-13  f: +2.662e-13  / ROOT___ 3

In [3]: mplot_s(f)
```
(matplotlib plot not shown!)

##### Square root of two
Here is a quick example of function inversion.
There is a choice of analysis (root finding) method:
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: newton_d(lambda x: x * x - 2.0, x0=1.0)
Out[1]: Result(method='NT', x=1.414213562373095, f=4.440892098500626e-16, δx=-1.570092458683775e-16, count=6, sense='_', mode='ROOT___')

In [2]: timeit(newton_d(lambda x: x * x - 2.0, x0=1.0))
25.4 µs ± 1.35 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)

In [3]: bisect_d(lambda x: x * x, xa=1.0, xb=2.0)
Out[3]: Result(method='BI', x=2.0, f=4.0, δx=0.0, count=101, sense='_', mode='ROOT___')

In [4]: timeit(bisect_d(lambda x: x * x, xa=1.0, xb=2.0))
342 µs ± 7.4 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```

##### Automatic Differentiation
Here we calculate _all_ the derivatives of a simple cubic in x, followed by its sensitivities to each parameter a, b, c.
```
$ ipython -i ipython.py
Python 3.6.8 (default, Jan 14 2019, 11:02:34)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
ad module loaded
plotters module loaded
ode_analysis module loaded

In [1]: a = Series.get(5, 3.0)

In [2]: b = Series.get(5, 5.0)

In [3]: c = Series.get(5, 7.0)

In [4]: x = Series.get(5, 2.0)

In [5]: print(a * x**3 - b * x**2 + c * x - 5)
+1.300e+01 +0.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [6]: print(a * x.var**3 - b * x.var**2 + c * x.var - 5)
+1.300e+01 +2.300e+01 +1.300e+01 +3.000e+00 +0.000e+00

In [7]: print(a.var * x**3 - b * x**2 + c * x - 5)
+1.300e+01 +8.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [8]: print(a * x**3 - b.var * x**2 + c * x - 5)
+1.300e+01 -4.000e+00 +0.000e+00 +0.000e+00 +0.000e+00

In [9]: print(a * x**3 - b * x**2 + c.var * x - 5)
+1.300e+01 +2.000e+00 +0.000e+00 +0.000e+00 +0.000e+00
```
