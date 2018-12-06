#!/usr/bin/env python3

# Example: ./series_test.py

from sys import stderr
from math import pi
from series import Series
from taylor import jet_c

print("", file=stderr)
sine, cosine = Series(jet_c(pi / 3.0, 7), diff=True).sin_cos
print(sine.derivatives, file=stderr)
print(cosine.derivatives, file=stderr)
tangent = Series(jet_c(pi / 4.0, 7), diff=True).tan
print(tangent.derivatives, file=stderr)
print("", file=stderr)

z = Series(jet_c(3, 7), diff=True)
print(z, file=stderr)
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
print(z + z, file=stderr)
print("", file=stderr)
print(z / 2, file=stderr)
print(2 / z, file=stderr)
print(z / z, file=stderr)
print("", file=stderr)
print(z ** 2, file=stderr)

print("", file=stderr)
print(Series(jet_c(2.0, 7), diff=True).sqrt, file=stderr)
