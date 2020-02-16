## Background

This project is mainly a collection of programs for evolving systems of ODEs using the Taylor Series Method (TSM), a rather old but poorly acknowledged technique based on forward mode Automatic Differentiation (AD).
TSM is a procedure for integrating ODEs using Taylor Series of arbitrary order, calculated to arbitrary precision (the former requires the latter in practice), using recurrence relations between time derivatives of increasing order.
It is thus (or should be!) a serious competitor to the fourth-order RK4 for the majority of practical cases.
For the uninitiated, here is a review of the method itself and it's history of repeated "discovery".

https://arxiv.org/abs/1111.7149

My work was inspired by the following paper and its associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

That software uses code generation to produce the recurrences and code to drive them, and it inspired me to try coding something by hand.
These programs are the result.
The main objective was to solve coupled nonlinear equations and investigate chaotic systems.
There is a maths.tex file containing a technical note describing the majority of the procedure (apart from the recurrence relations which are obtained from the source listed below).

The Taylor recurrence rules (the global "t-functions" in c and Python) generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order calculations.
The functions provided cover the basic algebraic operations (+ - * /), and also include several common functions:
* abs
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, a is a scalar)
* ln

The recurrence relations used here are derived in http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).
There are also factories for derivative "jets", and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, properly called,  are all that is needed to "integrate" systems of ODEs.
There is a fairly extensive collection of nonlinear ODEs already implemented, in the file tsm.py.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and others.

As part of the work to verify my implementation of these recurrence rules, I have added a demonstration of using Taylor series to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Additionally it is designed for finding inverse values (where real solutions exist) of complicated functions, not just their roots.

The "higher-level" ad_functions (or the Series class in Python) generate Taylor series "jets" in one go, so are only useful for univariate functions.
There is an overloaded operator "~" which extracts the actual derivative values from the taylor series coefficents, as well as additional operators for (negation and **).
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
There are three main areas of application for the code:
* solving nonlinear ODEs
* plotting functions and their (higher) derivatives, with solution, turning-point, and inflection analysis
* interactive investigations in the Python console

The playground.py script provides the analysis of a "model" function's derivatives along with the value itself, across a range of the input variable.
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

## Graph Plotting

There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python executables, which share output "formats").
In the example invocations given below, communication between the executable and plotting script uses a Unix pipe.
The dependencies are:
* matplotlib for 2D graphs and progressive graphs
* pi3d for 3D progressive trajectories (the visual python implementation is still distributed but is now considered legacy)
* gnuplot for static 3D trajectories

## Build Environment (Debian/Ubuntu/Raspbian)

OS level requirements:
```
sudo apt install build-essential mesa-utils-extra python3-dev libmpc-dev libatlas-base-dev virtualenvwrapper gnuplot-x11
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground
cd ODE-Playground
```
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
c Build (GCC or Clang)
```
$ ./build
$ ./build clang
```

## Running Python Tests
Testing Requirements
* pytest
* pytest_cov
* mutmut

Most of the code is covered several times over.
Tests have been mutation tested with mutmut.
```
pytest --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py -v
```

## Running c Tests
A successful test run creates no output:
```
$ ./build [clang]
<build output>
$ ./ad-test-dbg 7 2 1 >/tmp/ad-test.txt; diff --context=1 /tmp/ad-test.txt ad-test.txt
$
```

## Solving and Plotting ODEs
This use case only involves calling the "t-functions" in tsm.py.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm.py and tsm-*.c for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.
To find some example invocations:
```
grep Example *.py tsm-*.c
```
Matplotlib progressive graph plotting in Python (second parameter ignored as python floats are fixed double precision):
```
./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```
or in c (notice absence of the model parameter!):
```
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D static trajectory plotting (gnuplot)

Write to a data file
```
./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 >/tmp/data
```
then plot it
```
echo "splot '/tmp/data' with lines" | gnuplot -p
```

3D progressive trajectory plotting (pi3d)
```
./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
```

## tsm.py Parameter Reference
tsm.py comprises a long "if" statement containing a "zoo" of pre-programmed ODE systems.

Parameter | Meaning
----------|-----------
1 | model (ODE) name
2 | significant figures for display
3 | 1 + order
4,5 | step size, number of steps
6,7,8 | x0, y0, z0
9+ | ODE parameters

c is mostly the same, except the "model" parameter, which is part of the executable name, is omitted from the list, and the executable itself is named tsm-"model".

## Interactivity

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
