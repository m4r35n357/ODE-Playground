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

The Taylor recurrence rules generate "jets" of derivatives (the Taylor Series coefficients) term by term using previously calculated lower order derivatives.
The functions provided cover the basic algebraic operations (+ - * /), and also include several common functions:
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)_sec(h)2
* pwr
* ln

There are also factories for derivative "jets", and an implementation of Horner's method for summing the Taylor Series.
These "low-level" functions, properly called,  are all that is needed to "integrate" systems of ODEs.

As part of the work verifying my implementation of these recurrence rules, I added a demonstration of using Taylor series to find roots, extrema and inflection points in single variable nonlinear equations, much like the Matlab implementation described here.
http://www.neidinger.net/SIAMRev74362.pdf

The "higher-level" functions generate Taylor series "jets" in one go, so are only useful for univariate functions.
As well as adding operators for (negation and **), there are functions for:
* sqr
* sqrt
* exp
* sin(h)_cos(h)
* tan(h)
* ln
* abs

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
./tsm_lorenz.py 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50
```

3D ODE plotting
```
./tsm_lorenz.py 160 100 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
./tsm_lorenz.py 160 100 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotTrajectory.py 3 0 1 2
```
