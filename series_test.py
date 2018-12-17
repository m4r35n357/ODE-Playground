#!/usr/bin/env python3

# Example: ./series_test.py

from math import pi
from series import Series
from taylor import t_jet

a = 3.0
b = 4.0
sine, cosine = Series(t_jet(7, pi / a), diff=True).sin_cos
print("sin(pi / {})".format(a))
print(~ sine)
print("cos(pi / {})".format(a))
print(~ cosine)
print("tan(pi / {})".format(b))
print(~ Series(t_jet(7, pi / 4.0), diff=True).tan)

a = 3
b = 4
print("z = {}".format(a))
z = Series(t_jet(7, a), diff=True)
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
print(~ Series(t_jet(7, b), diff=True).sqrt)
