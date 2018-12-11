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
* pwr
* ln

The recurrence relations used here are derived in http://aimsciences.org/journals/displayPaperPro.jsp?paperID=9241 (open access).
There are also factories for derivative "jets", and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, properly called,  are all that is needed to "integrate" systems of ODEs.
There is a fairly extensive collection of nonlinear ODEs already implemented, either in c or Python, or both in most cases.

As part of the work verifying my implementation of these recurrence rules, I have added a demonstration of using Taylor series to implement Newton's method along the lines of the Matlab implementation described here http://www.neidinger.net/SIAMRev74362.pdf.
The demo can be used to find roots (and also extrema and inflection points by "extending" Newton to higher derivatives) in single variable nonlinear equations.


The "higher-level" functions generate Taylor series "jets" in one go, so are only useful for univariate functions.
As well as adding operators for (negation and **), there are functions for:
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)
* ln
* abs

Using these higher level functions, Newton's method is implemented trivially, but I have also provided implementations of bisection and arbitrary degree Householder methods for comparison.
The playground.py script generates the first 12 derivatives of a "model" function along with the value itself, across a range of the input variable.
Optionally it will analyse that range for roots, extrema and inflection points using the lower order derivatives.
Here is a seventh-degree polynomial model by way of example, and a trigonometric identity as a check of the sin and sqr functions:
```python
def septic(a, value):
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value

def trig(a, value):
    return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
```
Operators on or between jets and numbers are implemented by overloading, and the functions are implemented as properties of a jet.
These two approaches produce a fairly readable format for coding the models.
The value parameter allows simple calculation of function inverse values, as well as roots.

## Plotting
There are also some plotting and graphing utilities written in Python 3, which is required for plotting (the data can come from either c or Python, which share output "formats").
The dependencies are:
* matplotlib for graphs
* vpython *or* pi3d for 3D trajectories

## c version
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

To test the root finding built on it:
```
./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 8 8000 >/dev/null
```

Matplotlib progressive ODE plotting
```
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D ODE plotting
```
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotTrajectory.py 3 0 1 2
```


## Python version

There is also a (single precision) roughly equivalent Python 3 version of the programs with built-in models which has no external dependencies.
Dependencies of Python 3 programs:
* None

However, it is an easy task to convert the Python solver to use gmpy2 (MPFR) if desired.
To find Python example invocations:
```
grep Example *.py
```

To test the root, extremum and inflection finding:
```
./playground.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 8 8000 >/dev/null
```

Matplotlib progressive ODE plotting
```
./tsm-mp.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
./tsm-mp.py 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D ODE plotting
```
./tsm_lorenz.py 160 100 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
./tsm_lorenz.py 160 100 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotTrajectory.py 3 0 1 2
```
