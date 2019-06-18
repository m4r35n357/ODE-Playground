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

The Taylor recurrence rules generate "jets" of derivatives (the Taylor Series coefficients) term by term using previously calculated lower order derivatives.
The functions provided cover the basic algebraic operations (+ - * /), and also include several common functions:
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* asin (Python only)
* acos (Python only)
* atan (Python only)
* pwr
* ln

The recurrence relations used here are derived in http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).
There are also factories for derivative "jets", and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, properly called,  are all that is needed to "integrate" systems of ODEs.
There is a fairly extensive collection of nonlinear ODEs already implemented, either in c or Python, or both in most cases.
The list includes systems due to Lorenz, Rossler, Thomas, Bouali, Rabinovitch-Fabrikant, Sprott, and others.

As part of the work to verify my implementation of these recurrence rules, I have added a demonstration of using Taylor series to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The solver demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.
Additionally it is designed for finding inverse values (where real solutions exist) of complicated functions, not just their roots.


The "higher-level" functions generate Taylor series "jets" in one go, so are only useful for univariate functions.
As well as adding operators for (negation and **), there are functions for:
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)
* asin (Python only)
* acos (Python only)
* atan (Python only)
* ln
* abs

Using these higher level functions, Newton's method is implemented trivially, but I have also provided implementations of bisection and arbitrary degree Householder (c only) methods for comparison.
The playground.py script generates the first 12 derivatives of a "model" function along with the value itself, across a range of the input variable.
Optionally it will analyse that range for roots, extrema and inflection points using the lower order derivatives.
Here is a seventh-degree polynomial model by way of example, and a trigonometric identity as a check of the sin and sqr functions:
```python
def septic(a, value):
    return (a + 7) * (5 + a) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value

def trig(a, value):
    return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
```
Operators on or between jets and numbers are implemented by overloading, and the functions are implemented as properties of a jet.
These two approaches produce a fairly readable format for coding the models.
The value parameter allows simple calculation of function inverse values, as well as roots.

## Plotting
There are some plotting and graphing utilities written in Python 3, (the data itself can come from either c or Python, which share output "formats").
The dependencies are:
* matplotlib for 2D graphs
* pi3d for 3D trajectories (the visual python implementation is still distributed but is now considered legacy)

## Getting started with Python/gmpy2

There is a Python 3 version of the ODE solver programs with built-in models.
Dependencies of Python 3 programs:
* gmpy2 (based on MPFR)
* matplotlib
* pi3d

All future development will now be in Python3 with arbitrary precision, and the c/MPFR version will be kept as a reference implementation.

Build environment (Debian/Ubuntu)
```
sudo apt install build-essential mesa-utils-extra python3-dev virtualenvwrapper
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
To find Python ODE example invocations:
```
grep Example *.py
```
It is also worth looking for examples in the c files as below.
To test the root, extremum and inflection finding:
```
./models.py 2 -8 8 1001 7 1e-9 1e-9 | ./plotMany.py 8 10 >/dev/null
```

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

## Parameter reference
For Python, the first parameter is a string identifying the ODE, so numbers here refer only to the numeric parameters.

Parameter | Meaning
----------|-----------
1,2 | precision in decimal places, order
3,4 | step size, number of steps
5,6,7 | x0, y0, z0
8+ | ODE parameters


## c usage
Dependencies of c programs:
* MPFR 4 or later! (otherwise stick to the Python version)

Build them using the command:
```
./build
```

The built c programs are all called tsm-[something].
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


