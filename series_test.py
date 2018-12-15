#!/usr/bin/env python3

# Example: ./series_test.py

from sys import stderr
from math import pi
from series import Series
from taylor import t_jet

print("", file=stderr)
sine, cosine = Series(t_jet(7, pi / 3.0), diff=True).sin_cos
print(sine.derivatives, file=stderr)
print(cosine.derivatives, file=stderr)
tangent = Series(t_jet(7, pi / 4.0), diff=True).tan
print(tangent.derivatives, file=stderr)
print("", file=stderr)

z = Series(t_jet(7, 3), diff=True)
print(z, file=stderr)
print(- z, file=stderr)
print("", file=stderr)
print(z + 2, file=stderr)
print(2 + z, file=stderr)
print(z + z, file=stderr)
print("", file=stderr)
print(z - 2, file=stderr)
print(2 - z, file=stderr)
print(z - z, file=stderr)
print("", file=stderr)
print(z * 2, file=stderr)
print(2 * z, file=stderr)
print(z * z, file=stderr)
print("", file=stderr)
print(z / 2, file=stderr)
print(2 / z, file=stderr)
print(z / z, file=stderr)
print("", file=stderr)
print(z ** 2, file=stderr)

print("", file=stderr)
print(Series(t_jet(7, 2.0), diff=True).sqrt, file=stderr)
