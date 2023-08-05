#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

from ctypes import *
from math import acos, factorial

class DualNumber(Structure):
    _fields_ = [("val", c_longdouble), ("dot", c_longdouble)]


dual = CDLL('./libdual.so')
dual.d_dual.restype = DualNumber
dual.d_var.restype = DualNumber

a = dual.d_dual(c_longdouble(7.0))
print(a.val, a.dot)

b = dual.d_var(c_longdouble(5.0))
print(b.val, b.dot)

class Pair(Structure):
    _fields_ = [("a", c_longdouble), ("b", c_longdouble)]


TaylorSeries = c_longdouble * 4

taylor = CDLL('./librecurrences.so')
taylor.t_mul.restype = c_longdouble
taylor.t_exp.restype = c_longdouble
taylor.t_sin_cos.restype = Pair
taylor.t_ln.restype = c_longdouble

a = TaylorSeries()
a[0] = c_longdouble(7.0)
a[1] = c_longdouble(1.0)
print(a)
for k in range(len(a)):
    print(a[k])

b = TaylorSeries()
b[0] = c_longdouble(5.0)
b[1] = c_longdouble(1.0)
print(b)
for k in range(len(b)):
    print(b[k])

c = TaylorSeries()
print(c)
for k in range(len(c)):
    c[k] = taylor.t_mul(a, b, k)
    print(c[k])

a = TaylorSeries()
a[0] = c_longdouble(1.0)
a[1] = c_longdouble(1.0)
c = TaylorSeries()
print(c)
for k in range(len(a)):
    taylor.t_exp(c, a, k)
    print(c[k] * factorial(k))

a = TaylorSeries()
a[0] = c_longdouble(1.0)
a[1] = c_longdouble(1.0)
c = TaylorSeries()
print(c)
for k in range(len(a)):
    taylor.t_ln(c, a, k)
    print(c[k] * factorial(k))

a = TaylorSeries()
a[0] = c_longdouble(acos(-1.0))
a[1] = c_longdouble(1.0)
c = TaylorSeries()
d = TaylorSeries()
print(c, d)

for k in range(len(a)):
    taylor.t_sin_cos(byref(c), byref(d), a, k)
    fac = factorial(k)
    print(c[k] * fac, d[k] * fac)

a = TaylorSeries()
a[0] = c_longdouble(acos(-1.0))
a[1] = c_longdouble(1.0)
c = TaylorSeries()
d = TaylorSeries()
print(c, d)
#e = Pair()

for k in range(len(a)):
    e = taylor.t_sin_cos(c, d, a, k)
    fac = factorial(k)
    print(e.a * fac, e.b * fac)
