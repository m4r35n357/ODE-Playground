
This is a collection of programs for evolving systems of ODE using the Taylor Series Method.  It was inspired by this paper and associated software:

https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

That software uses code generation, and it inspired me to try coding something by hand.  These programs are the result.  The main objective was to solve coupled nonlinear equations and investigate chaotic systems.  As part of the verification I have added a demonstration of ustin Taylor series to find roots in single variable nonlinear equations.

Dependencies of c programs:
MPFR 4 or later! (otherwise stick to the Python version)

Build them using the command:

./build

The built programs are all called tsm-<NAME>... .  Each source file should contain an example invocation near the top.

There is also a (single precision) roughly equivalent Python 3 version of the solver with built-in models which has no external dependencies.  There are also some plotting and graphing utilities written in Python 3 which need:

vpython
pi3d
matplotlib

It is an easy task to convert the Python solver to use gmpy2 (MPFR)

To test the top level Taylor series operation:
./ad-test-dbg 7 2 1

To test the root finding built on it:
./ad-test-newton-dbg 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 8 8000 >/dev/null

To generate and plot an ODE using lower level operations:

(Python version with matplotlib progressive plotting)
./tsm-mp.py lorenz 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -30 50

(c version with Pi3D trajectory plot)
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
(c version with visula python trajectory plot)
./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3 | ./plotTrajectory.py 3 0 1 2

For a few more ideas:

grep Example *.c

