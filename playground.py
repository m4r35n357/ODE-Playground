#!/usr/bin/env python3

# Example: ./playground.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null

from sys import argv, stderr
from math import pi
from series import Series, newton, householder
from taylor import jet_0, jet_c


def cosx_x3(a, value):
    return a.cos - a * a * a - value


def septic(a, value):
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def playground(a, value):
    # return (a.sqr)**0.5 - value
    return a.sqrt - a ** 0.5 - value


N_MAX = 13
n = int(argv[1])
x0 = float(argv[2])
x1 = float(argv[3])
steps = int(argv[4])
target = float(argv[5])
fun = playground

x_step = (x1 - x0) / steps
w_x = Series(jet_c(x0, N_MAX), diff=True)
w_f = Series(jet_0(N_MAX))
for k in range(steps):
    w_x.jet[0] = x0 + k * x_step
    w_f = fun(w_x, 0.0)
    print("{:.6e} {}".format(w_x.jet[0], w_f.derivatives))
    if k > 0:
        # noinspection PyUnboundLocalVariable
        if f_prev * w_f.jet[0] < 0.0:
            print("Bracketed root, solving", file=stderr)
            if n == 2:
                print("using Newton's method", file=stderr)
                newton(fun, w_x.jet[0])
            else:
                print("using Householder's method", file=stderr)
                householder(fun, w_x.jet[0], n)
    x_prev = w_x.jet[0]
    f_prev = w_f.jet[0]


print("", file=stderr)
sine, cosine = Series(jet_c(pi / 3.0, 7), diff=True).sin_cos
print(sine, file=stderr)
print(sine.derivatives, file=stderr)
print(cosine, file=stderr)
print(cosine.derivatives, file=stderr)
print("", file=stderr)

tangent, secant2 = Series(jet_c(pi / 4.0, 7), diff=True).tan_sec2
print(tangent, file=stderr)
print(tangent.derivatives, file=stderr)
print(secant2, file=stderr)
print(secant2.derivatives, file=stderr)
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
