#!/usr/bin/env python3

# Example: ./playground.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null

from sys import argv, stderr
from series import Series, newton, householder
from taylor import jet_0, jet_c


def cosx_x3(a, value):
    return a.cos - a * a * a - value


def septic(a, value):
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def playground(a, value):
    return ((2 * a ** 2) / 3) * (3 / (2 * a ** 2)) - value


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
