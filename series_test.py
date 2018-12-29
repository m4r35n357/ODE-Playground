#!/usr/bin/env python3

# Example: ./series_test.py

from math import pi
from gmpy2 import get_context
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from series import Series
from taylor import t_jet, to_mpfr

# noinspection PyArgumentList
a = to_mpfr("3.0")
# noinspection PyArgumentList
b = to_mpfr("4.0")
sine, cosine = Series(t_jet(7, pi / a), variable=True).sin_cos
print("sin(pi / {})".format(a))
print(~ sine)
print("cos(pi / {})".format(a))
print(~ cosine)
print("tan(pi / {})".format(b))
print(~ Series(t_jet(7, pi / 4.0), variable=True).tan)

# noinspection PyArgumentList
a = to_mpfr(3)
# noinspection PyArgumentList
b = to_mpfr(4)
print("z = {}".format(a))
z = Series(t_jet(7, a), variable=True)
print("+, -")
print(+ z)
print(- z)
print("+")
print(z + b)
print(b + z)
print(z + z)
print("-")
print(z - b)
print(b - z)
print(z - z)
print("*")
print(z * b)
print(b * z)
print(z * z)
print("/")
print(z / b)
print(b / z)
print(z / z)
print("**b")
print(z ** b)

print("sqrt")
print(~ Series(t_jet(7, b), variable=True).sqrt)
