
This is a collection of programs for evolving systems of ODE using the Taylor Series Method.
It was inspired by this paper and associated software:
https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

That software uses code generation, and it inspired me to try coding something by hand.
These programs are the result.
The main objective was to solve coupled nonlinear equations and investigate chaotic systems.

As part of the work verifying the lower layers, I added a demonstration of using Taylor series to find roots, extrema and inflection points in single variable nonlinear equations, like here.
http://www.neidinger.net/SIAMRev74362.pdf

# Plotting
There are also some plotting and graphing utilities written in Python 3, which is required for plotting (the data can come from either c or Python, which share output "formats").
The dependencies are:
* matplotlib for graphs
* vpython *or* pi3d for 3D trajectories

# c version
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


# Python version

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
