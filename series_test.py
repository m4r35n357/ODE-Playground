#!/usr/bin/env python3

# Example: ./series_test.py

from math import pi
from gmpy2 import get_context
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from series import Series
from taylor import t_jet, to_mpfr

a = to_mpfr("3.0")
b = to_mpfr("4.0")
sine, cosine = Series(t_jet(7, pi / a), variable=True).sin_cos
print("sin(pi / {})".format(a))
print(~ sine)
print("cos(pi / {})".format(a))
print(~ cosine)
print("tan(pi / {})".format(b))
print(~ Series(t_jet(7, pi / 4.0), variable=True).tan)
c = to_mpfr("0.5")
arc_sine = Series(t_jet(7, c), variable=True).asin
print("asin(0.5)".format(c))
print(arc_sine)
arc_cosine = Series(t_jet(7, c), variable=True).acos
print("acos(0.5)".format(c))
print(arc_cosine)
arc_tangent = Series(t_jet(7, c), variable=True).atan
print("atan(0.5)".format(c))
print(arc_tangent)

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
