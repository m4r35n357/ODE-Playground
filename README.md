## Background

This project is mainly a collection of programs in c and Python for evolving systems of ODEs using the Taylor Series Method (TSM), a rather old but poorly acknowledged technique based on forward mode Automatic Differentiation (AD).
TSM is a procedure for integrating ODEs using Taylor Series of arbitrary order, calculated to arbitrary precision (the former requires the latter in practice), using recurrence relations between time derivatives of increasing order.
It is therefore (or should be, if it were better known!) a serious competitor to the fourth-order RK4 for the vast majority of cases.

The code itself is tiny (as are the dynamic executables) and has been partly developed on, and is suitable for running on, a Raspberry Pi 4 computer.
The c code supports arbitrary precision, and is most suited to solving systems of ODEs to high accuracy.
The Python code uses float precision, and is most suited to interactive analysis and plotting of functions and their derivatives.

For the uninitiated, here is a review of the TSM itself and its history of repeated "discovery" and re-branding.
https://arxiv.org/abs/1111.7149

My work was inspired by parts of the following paper and its associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

That software uses code generation to produce the recurrences and code to drive them, and it inspired me to try coding something more direct by hand.
These programs are the result.
My primary aim was to be able to solve coupled nonlinear equations and investigate chaotic systems, without relying on "black-box" ODE solvers.
The header file taylor-ode.h contains a terse but complete description of the Taylor Series Method and the derivations of the recurrences that enable analysis of complex systems.

The recurrence rules (the "t-functions" in c and Python) are the key to calculating high order derivatives accurately, without needing finite differences.
They generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order values.
The functions provided cover the basic algebraic operations on Taylor Series (+ - * /), and also include several common functions:
* abs
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, a is a scalar)
* ln

The recurrence relations used here are derived along the lines of (amongst other sources) http://www2.math.uni-wuppertal.de/wrswt/preprints/prep_05_4.pdf and http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).

There are also convenient factories for generating derivative "jets" of arbitrary order, and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, when properly called,  are all that is needed to solve systems of ODEs.
There is a fairly extensive collection of nonlinear ODE examples already implemented, in the file tsm.py, and in the tsm-\*-*.c files.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and many others.

As part of the work to verify my own implementation of these recurrence rules, I have added a demonstration of using Taylor series to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Additionally it is designed for finding inverse values (where real solutions exist) of complicated functions, not just their roots.

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
* pwr (f(x)^a, a is a scalar)
* ln

Using these higher level functions, Newton's method is implemented trivially, but I have also provided an implementations of the bisection method for comparison.
In summary, there are three main areas of application for the code:
* solving nonlinear ODEs
* plotting functions and their (higher) derivatives, with solution, turning-point, and inflection analysis
* interactive investigations in the Python console

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
pytest --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py -v
```

#### c Build (GCC or Clang)
```
$ ./build
$ ./build clang
```

#### Running c Tests
##### Newer tests
```
$ ./libad-test-dbg 9 32 20 2 1e-18 1
```
The final parameter can be set to 0 (or left absent) for a summary, or 2 for full detail.

Parameter | Meaning
----------|-----------
1 | (approximate) precision in decimal places
2 | order of Taylor Series
3 | x value
4 | tolerance
5 | debug level (optional)

##### Older tests
```
$ ./ad-test-dbg 7 2 1
```
Output originally designed for visual checking, or see below how to do it automatically (hard-coded to quad precision).

Parameter | Meaning
----------|-----------
1 | order of Taylor Series
2 | x value
3 | y value

##### Big build and test command:
```
time -p ./build && ./ad-test-dbg 7 2 1 >/tmp/ad-test.txt; diff --context=1 /tmp/ad-test.txt ad-test.txt && ./libad-test-dbg 9 32 20 2 1e-18 && echo OK
```

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
$ ./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```
or in c (notice absence of the model parameter!):
```
$ ./tsm-lorenz-dbg 15 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -30 50
```

#### tsm.py and tsm-\*-\* Parameter Reference
The Python ODE solver, tsm.py, comprises a long "if" statement containing a "zoo" of pre-programmed ODE systems.
tsm-\*-\* represents the individual c solvers.

##### Python (float precision)
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
$ ./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 >/tmp/data
```
then plot it
```
$ gnuplot -p -e "splot '/tmp/data' with lines"
```

#### 3D progressive trajectory plotting (pi3d)
```
$ ./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plot3d.py
```

## Clean Numerical Simulation (CNS)

As a rough guide to the accuracy of a solution, it can be compared with a "better" solution (one made with "better" solver parameters), and discarded at the point where the solutions diverge.
The two simulations are run in parallel processes, but obviously the "better" solution takes longer.

300 time units
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

## Sensitivity to variation in initial conditions

The simulation is run six times in parallel processes, with each x, y, z initial condition perturbed by +/- the parameter to ic.
```
$ ./ic .001 ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
x+ ./tsm-lorenz-dbg 16 10 .01 10001 -15.799 -17.48 35.64 10 28 8 3
x- ./tsm-lorenz-dbg 16 10 .01 10001 -15.801 -17.48 35.64 10 28 8 3
y+ ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.479 35.64 10 28 8 3
y- ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.481 35.64 10 28 8 3
z+ ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.641 10 28 8 3
z- ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.639 10 28 8 3
```
(3D plot not shown)

## Interactivity with Python
This use case involves calling the Series methods and operators from ad.py, and the functions from plotters.py, from a Python interpreter.

##### Higher level function analysis
Plotting function and derivatives, together with root and turning point analysis:
```
$ ipython3
Python 3.7.1 (default, Oct 22 2018, 11:21:55)
Type 'copyright', 'credits' or 'license' for more information
IPython 7.2.0 -- An enhanced Interactive Python. Type '?' for help.
PyDev console: using IPython 7.2.0
Python 3.7.1 (default, Oct 22 2018, 11:21:55)
[GCC 8.2.0] on linux
from plotters import *
Backend TkAgg is interactive backend. Turning interactive mode on.
ad module loaded
playground module loaded
plotters module loaded
f = lambda a: (a.exp + (a.sqr - 4.0).exp).ln
scan(f)
BI  x: -1.962e+00  δx: -2.387e-10  f: +4.026e-10  \ ROOT___ 26
BI  x: -1.312e+00  δx: -4.773e-10  f: +2.023e-10  / MIN_MAX 25
BI  x: -1.849e-02  δx: -9.546e-10  f: -2.642e-10  / ROOT___ 24
mplot(f)
```
(matplotlib plot not shown!)

##### Square root of two
Here is a quick example of function inversion.
There is a choice of analysis (root finding) method:
```
$ ipython3
...
In [1]: from ad import *
   ...: from playground import *
ad module loaded
playground module loaded

In [2]: bisect(lambda x: x**2 - 2, xa=0.0, xb=2.0)
Out[2]: Result(method='BI', x=1.414213562372879, f=-6.108447081487611e-13, δx=4.547473508864641e-13, count=42, sense='_', mode='ROOT___')

In [3]: newton(lambda x: x**2 - 2, x0=1.0)
Out[3]: Result(method='NT', x=1.414213562373095, f=4.440892098500626e-16, δx=-1.570092458683775e-16, count=6, sense='_', mode='ROOT___')

In [4]: timeit(bisect(lambda x: x**2 - 2, xa=0.0, xb=2.0))
600 µs ± 10.4 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)

In [5]: timeit(newton(lambda x: x**2 - 2, x0=1.0))
76 µs ± 1.63 µs per loop (mean ± std. dev. of 7 runs, 10000 loops each)
```

##### Automatic Differentiation
Here we calculate _all_ the derivatives of a simple cubic in x, followed by its sensitivities to each parameter a, b, c.
```
$ ipython3
...
In [3]: from ad import *
ad module loaded

In [4]: a = Series.get(5, 3.0)

In [5]: b = Series.get(5, 5.0)

In [6]: c = Series.get(5, 7.0)

In [7]: x = Series.get(5, 2.0)

In [8]: print(a * x**3 - b * x**2 + c * x - 5)
+1.300000e+01 +0.000000e+00 +0.000000e+00 +0.000000e+00 +0.000000e+00

In [9]: print(a * x.var**3 - b * x.var**2 + c * x.var - 5)
+1.300000e+01 +2.300000e+01 +1.300000e+01 +3.000000e+00 +0.000000e+00

In [10]: print(a.var * x**3 - b * x**2 + c * x - 5)
+1.300000e+01 +8.000000e+00 +0.000000e+00 +0.000000e+00 +0.000000e+00

In [11]: print(a * x**3 - b.var * x**2 + c * x - 5)
+1.300000e+01 -4.000000e+00 +0.000000e+00 +0.000000e+00 +0.000000e+00

In [12]: print(a * x**3 - b * x**2 + c.var * x - 5)
+1.300000e+01 +2.000000e+00 +0.000000e+00 +0.000000e+00 +0.000000e+00

```
