
This is a collection of programs for evolving systems of ODE using the Taylor Series Method.  It was inspired by this paper and associated software:

https://web.ma.utexas.edu/users/mzou/taylor/taylor.pdf

That software uses code generation, and it inspired me to try coding something by hand.  These programs are the result.

Dependencies:
MPFR

Build them using the command:

./build

The built programs are all called tsm-<NAME>... .  Each source file should contain an example invocation near the top.

There is also a (single precision) roughly equivalent Python version which has no external dependencies.

There are also some plotting and graphing utilities written in Python 3 which need:

vpython
pi3d
matplotlib

