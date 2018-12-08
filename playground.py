#!/usr/bin/env python3

# Example: ./playground.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null

from sys import argv, stderr
from series import SolveMode, Series, bisect, newton, householder
from taylor import jet_0, jet_c


def cosx_x3(a, value):
    return a.cos - a * a * a - value


def septic(a, value):
    return (a + 7) * (a + 5) * (a + 2) * a * (a - 1) * (a - 3) * (a - 6) - value


def playground(a, value):
    # return (2 * a).sin - 2 * a.sin * a.cos - value
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr - value
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
    # return a.sqr.sqrt - value
    return a.exp.ln - value


N_MAX = 13
n = int(argv[1])
x0 = float(argv[2])
x1 = float(argv[3])
steps = int(argv[4])
target = float(argv[5])

fun = septic

x_step = (x1 - x0) / steps
w_x = Series(jet_c(x0, N_MAX), diff=True)
w_f = Series(jet_0(N_MAX))
if n != 0:
    if n == 1:
        print("Bisection", file=stderr)
    elif n == 2:
        print("Newton's method", file=stderr)
    else:
        print("Householder's method of degree {}".format(n - 1), file=stderr)
for k in range(steps):
    w_x.jet[0] = x0 + k * x_step
    w_f = fun(w_x, target)
    print("{:.6e} {}".format(w_x.jet[0], w_f.derivatives))
    if n != 0:
        if k > 0:
            # noinspection PyUnboundLocalVariable
            if f_prev * w_f.jet[0] < 0.0:
                print("ROOT", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    bisect(fun, w_x.jet[0], x_prev, target=target)
                elif n == 2:
                    newton(fun, w_x.jet[0], target=target)
                else:
                    householder(fun, w_x.jet[0], n, target=target)
            # noinspection PyUnboundLocalVariable
            if f_dash_prev * w_f.jet[1] < 0.0:
                if f_dash_prev > w_f.jet[1]:
                    print("MAXIMUM", file=stderr)
                else:
                    print("MINIMUM", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    bisect(fun, w_x.jet[0], x_prev, mode=SolveMode.EXTREMUM)
                elif n == 2:
                    newton(fun, w_x.jet[0], mode=SolveMode.EXTREMUM)
                else:
                    householder(fun, w_x.jet[0], n, mode=SolveMode.EXTREMUM)
            # noinspection PyUnboundLocalVariable
            if f_dash_dash_prev * w_f.jet[2] < 0.0:
                print("INFLECTION", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    bisect(fun, w_x.jet[0], x_prev, mode=SolveMode.INFLECTION)
                elif n == 2:
                    newton(fun, w_x.jet[0], mode=SolveMode.INFLECTION)
                else:
                    householder(fun, w_x.jet[0], n, mode=SolveMode.INFLECTION)
    x_prev = w_x.jet[0]
    f_prev = w_f.jet[0]
    f_dash_prev = w_f.jet[1]
    f_dash_dash_prev = w_f.jet[2]
