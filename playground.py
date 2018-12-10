#!/usr/bin/env python3

#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

# Example: ./playground.py 2 -8 8 1001 0 1e-9 1e-9 | ./plotMany.py 8 50000 >/dev/null

from sys import argv, stderr
from series import Solver, Series, bisect, newton, householder
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


def print_output(output):
    print("    {:3d} {:22.15e} {:10.3e} {:10.3e}".format(*output), file=stderr)


N_MAX = 13
n = int(argv[1])
x0 = float(argv[2])
x1 = float(argv[3])
steps = int(argv[4])
target = float(argv[5])
f_tol = float(argv[6])
x_tol = float(argv[7])

fun = septic

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
    w_x.jet[0] = x0 + k * (x1 - x0) / steps
    w_f = fun(w_x, target).derivatives
    print("{:.6e} {}".format(w_x.jet[0], w_f))
    if n != 0:
        if k > 0:
            # noinspection PyUnboundLocalVariable
            if f_prev * w_f.jet[0] < 0.0:
                print("ROOT", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    print_output(bisect(fun, w_x.jet[0], x_prev, f_tol, x_tol, target=target))
                elif n == 2:
                    print_output(newton(fun, w_x.jet[0], f_tol, x_tol, target=target))
                else:
                    print_output(householder(fun, w_x.jet[0], n, f_tol, x_tol, target=target))
            # noinspection PyUnboundLocalVariable
            if f_dash_prev * w_f.jet[1] < 0.0:
                if f_dash_prev > w_f.jet[1]:
                    print("MAXIMUM", file=stderr)
                else:
                    print("MINIMUM", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    print_output(bisect(fun, w_x.jet[0], x_prev, f_tol, x_tol, mode=Solver.EXTREMUM))
                elif n == 2:
                    print_output(newton(fun, w_x.jet[0], f_tol, x_tol, mode=Solver.EXTREMUM))
                else:
                    print_output(householder(fun, w_x.jet[0], n, f_tol, x_tol, mode=Solver.EXTREMUM))
            # noinspection PyUnboundLocalVariable
            if f_dash_dash_prev * w_f.jet[2] < 0.0:
                if f_dash_dash_prev > w_f.jet[2]:
                    print("+INFLECTION", file=stderr)
                else:
                    print("-INFLECTION", file=stderr)
                if n == 1:
                    # noinspection PyUnboundLocalVariable
                    print_output(bisect(fun, w_x.jet[0], x_prev, f_tol, x_tol, mode=Solver.INFLECTION))
                elif n == 2:
                    print_output(newton(fun, w_x.jet[0], f_tol, x_tol, mode=Solver.INFLECTION))
                else:
                    print_output(householder(fun, w_x.jet[0], n, f_tol, x_tol, mode=Solver.INFLECTION))
    x_prev = w_x.jet[0]
    f_prev = w_f.jet[0]
    f_dash_prev = w_f.jet[1]
    f_dash_dash_prev = w_f.jet[2]
