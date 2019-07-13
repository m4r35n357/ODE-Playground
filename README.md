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

The Taylor recurrence rules (the global "t-functions") generate "jets" of Taylor Series coefficients iteratively, term by term, using previously calculated lower order calculations.
The functions provided cover the basic algebraic operations (+ - * /), and also include several common functions:
* abs
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr (f(x)^a, a is a scalar)
* ln

The recurrence relations used here are derived in http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).
There are also factories for derivative "jets", and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, properly called,  are all that is needed to "integrate" systems of ODEs.
There is a fairly extensive collection of nonlinear ODEs already implemented, either in c or Python, or both in most cases.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and others.

As part of the work to verify my implementation of these recurrence rules, I have added a demonstration of using Taylor series to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Additionally it is designed for finding inverse values (where real solutions exist) of complicated functions, not just their roots.


The "higher-level" functions (in the Series and Dual classes) generate Taylor series "jets" in one go, so are only useful for univariate functions.
There is an overloaded operator "~" which extracts the actual derivative values from the taylor series coefficents, as well as additional operators for (negation and **).
The \*\* (power) operator caters for f(x)^a, a^f(x) and f1(x)^f2(x), subject to domain limitations on f(x) and the scalar a.
There are also functions for:
* abs
* exp
* sin(h)_cos(h) (Python gmpy2 or Series class only)
* sin(h)
* cos(h)
* tan(h)
* ln

Using these higher level functions, Newton's method is implemented trivially, but I have also provided implementations of bisection and arbitrary degree Householder (c only) methods for comparison.
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
    return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr
```
Operators on or between jets and numbers are implemented by overloading, and the functions are implemented as _properties_ of a jet.
These two approaches produce a fairly readable format for coding the models.
The value parameter allows simple calculation of function inverse values, as well as roots.

## Graph Plotting
There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python executables, which share output "formats").
In the example invocations given below, communication between the executable and plotting script uses a Unix pipe.
The dependencies are:
* matplotlib for 2D graphs
* pi3d for 3D trajectories (the visual python implementation is still distributed but is now considered legacy)

## Getting started with Python/gmpy2

There is a Python 3 version of the ODE solver programs with built-in models.
Dependencies of Python 3 programs:
* gmpy2 (based on MPFR)
* matplotlib
* pi3d

Main future development will now be in Python3 with arbitrary precision (although currently there is a float-only branch), and the c/MPFR version will be kept as a reference implementation.

Build environment (Debian/Ubuntu/Raspbian)
```
sudo apt install build-essential mesa-utils-extra python3-dev libmpc-dev libatlas-base-dev virtualenvwrapper
```
This is a typical virtual environment setup:
```
mkvirtualenv --python /usr/bin/python3 ad
pip install gmpy2 matplotlib pillow pi3d
```
Download:
```
git clone https://github.com/m4r35n357/ODE-Playground
cd ODE-Playground
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

## Solving ODEs
This use case only involves calling the "t-functions" in tsm.py.
No differentiation happens in these functions (they only implement the recurrence relations); it is the responsibility of the calling program to organize this properly.
Refer to tsm.py for a varied selection of examples, including several from https://chaoticatmospheres.com/mathrules-strange-attractors.
To find Python ODE example invocations:
```
grep Example *.py
```
It is also worth looking for examples in the c files as below.

Matplotlib progressive graph plotting
```
./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D ODE plotting
```
./tsm.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
```

Double pendulum (see code for parameters)
```
./double.py 8 0.1 1000 1 1 1 1 .1 0 .1 0 | ./plotPi2d.py
```
## tsm.py Parameter reference
tsm.py comprises a long "if" statement containing a "zoo" of pre-programmed ODE systems.

Parameter | Meaning
----------|-----------
1 | model (ODE) name
2,3 | precision in decimal places, order
4,5 | step size, number of steps
6,7,8 | x0, y0, z0
9+ | ODE parameters


## Analysing functions
To test the function and multi-derivative plotting (Dual and Series representations), together with root, extremum and inflection finding:
```
./series_test.py 7 -8 8 1001 | ./plotMany.py 8 10 >/dev/null
```
```
$ ./models.py 1 -8 8 1001 13 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
ad module loaded
functions module loaded
playground module loaded
Bisection
ResultType(count=36, sense='-', mode='ROOT', x=mpfr('-1.961750534675549715757369995117187499999999999999999999999999999999999987',236), f=mpfr('-5.805043894073777177930058241509277124986288774641499850942009520578458612e-13',236), dx=mpfr('-4.656612873077392578124999999999999999999999999999999999999997102182694775e-13',236))
ResultType(count=2, sense='+', mode='EXTREMUM', x=mpfr('-1.320000000000000000000000000000000000000000000000000000000000000000000004',236), f=mpfr('-0.9895699319474545589789920831749399410839476209598406981422940090028371208',236), dx=mpfr('0.0',236))
ResultType(count=35, sense='+', mode='ROOT', x=mpfr('-0.01849182777944952249526977539062499999999999999999999999999999999999997151',236), f=mpfr('-9.131427064507876845535703720933450904719073854466606849410430432029352248e-13',236), dx=mpfr('-9.313225746154785156249999999999999999999999999999999999999999864164813818e-13',236))
```


## Interactivity
Here is a quick example of interactive use; function inversion.
```
$ ipython3 
Python 3.7.1 (default, Oct 22 2018, 11:21:55) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.2.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: from ad import *                                                        
ad module loaded

In [2]: from playground import *                                                
playground module loaded

In [3]: newton(lambda x: (x**2), x0=1.0, target=2.0)                            
Out[3]: ResultType(count=7, sense='', mode='ROOT', x=1.414213562373095, f=2.0000000000000004, dx=-1.5700924586837747e-16)

In [4]: bisect(lambda x: (x**2), xa=1.4, xb=1.5, target=2.0)                    
Out[4]: ResultType(count=38, sense='', mode='ROOT', x=1.4142135623733338, f=2.0000000000006755, dx=7.276401703393276e-13)
```
Here we differentiate a simple function of three variables with respect to each one.
This is done twice; first using the Dual class and then using the Series class.
The latter provides derivatives of higher order than the first.
```
$ ipython3
...
In [1]: from ad import *
ad module loaded

In [2]: a = Dual.get(3.0)

In [3]: b = Dual.get(5.0)

In [4]: c = Dual.get(7.0)

In [5]: print(a * c**2 - b * c)
+1.120000e+02 +0.000000e+00

In [6]: print(a.var * c**2 - b * c)
+1.120000e+02 +4.900000e+01

In [7]: print(a * c**2 - b.var * c)
+1.120000e+02 -7.000000e+00

In [8]: print(a * c.var**2 - b * c.var)
+1.120000e+02 +3.700000e+01

In [9]: a = Series.get(3, 3.0)

In [10]: b = Series.get(3, 5.0)

In [11]: c = Series.get(3, 7.0)

In [12]: print(a * c**2 - b * c)
+1.120000e+02 +0.000000e+00 +0.000000e+00 

In [13]: print(a.var * c**2 - b * c)
+1.120000e+02 +4.900000e+01 +0.000000e+00 

In [14]: print(a * c**2 - b.var * c)
+1.120000e+02 -7.000000e+00 +0.000000e+00 

In [15]: print(a * c.var**2 - b * c.var)
+1.120000e+02 +3.700000e+01 +3.000000e+00 
```
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

## c usage

This documentation is becoming obsolete, but the programs themselves will be retained.
Dependencies of c programs:
* MPFR 4 or later! (otherwise stick to the Python version)

Build them using the command:
```
./build
```

The built c programs are all called tsm-[model].
Each source file should contain an example invocation near the top.

To see them all:

```
grep Example *.c
```
To test the top level Taylor series operation:
```
./ad-test-dbg 7 2 1
```

To test the root finding:
```
./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
```

Matplotlib progressive graph plotting
```
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D ODE plotting
```
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
```
## tsm-[model] Parameter reference
Parameter | Meaning
----------|-----------
1,2 | precision in decimal places, order
3,4 | step size, number of steps
5,6,7 | x0, y0, z0
8+ | ODE parameters


