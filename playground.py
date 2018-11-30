#!/usr/bin/env python3

from math import pi
from series import newton, householder, Series
from taylor import jet_c


def fun(a, value):
    return a.sin_cos[1] - a * a * a - value


newton(fun, 0.5, target=0.0)
print("")

householder(fun, 0.5, 4, target=0.0)
print("")

sine, cosine = Series(jet_c(pi / 3.0, 7), diff=True).sin_cos
print(sine)
print(sine.derivatives)
print(cosine)
print(cosine.derivatives)
print("")

tangent, secant2 = Series(jet_c(pi / 4.0, 7), diff=True).tan_sec2
print(tangent)
print(tangent.derivatives)
print(secant2)
print(secant2.derivatives)
print("")

z = Series(jet_c(3, 7), diff=True)
print(z)
print("")
print(z + 2)
print(2 + z)
print(z + z)
print("")
print(z - 2)
print(2 - z)
print(z - z)
print("")
print(z * 2)
print(2 * z)
print(z + z)
print("")
print(z / 2)
print(2 / z)
print(z / z)
print("")
print(z ** 2)
