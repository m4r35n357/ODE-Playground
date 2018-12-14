#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv, stderr
from taylor import jet_c
from series import Solver, Series, bisect, newton, householder


def _print_output(output):
    print("    {:3d} {:22.15e} {:10.3e} {:10.3e}".format(*output), file=stderr)


def analyze(model):
    n_max = 13
    n = int(argv[1])
    x0 = float(argv[2])
    x1 = float(argv[3])
    steps = int(argv[4])
    target = float(argv[5])
    f_tol = float(argv[6])
    x_tol = float(argv[7])
    w_x = Series(jet_c(x0, n_max), diff=True)
    if n != 0:
        if n == 1:
            print("Bisection", file=stderr)
        elif n == 2:
            print("Newton's method", file=stderr)
        else:
            print("Householder's method of degree {}".format(n - 1), file=stderr)
    for k in range(steps):
        w_x.jet[0] = x0 + k * (x1 - x0) / steps
        w_f = model(w_x, target).derivatives
        print("{:.6e} {}".format(w_x.jet[0], w_f))
        if n != 0:
            if k > 0:
                # noinspection PyUnboundLocalVariable
                if f_prev * w_f.jet[0] < 0.0:
                    print("ROOT", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        _print_output(bisect(model, w_x.jet[0], x_prev, f_tol, x_tol, target=target))
                    elif n == 2:
                        _print_output(newton(model, w_x.jet[0], f_tol, x_tol, target=target))
                    else:
                        _print_output(householder(model, w_x.jet[0], n, f_tol, x_tol, target=target))
                # noinspection PyUnboundLocalVariable
                if f_dash_prev * w_f.jet[1] < 0.0:
                    if f_dash_prev > w_f.jet[1]:
                        print("MAXIMUM", file=stderr)
                    else:
                        print("MINIMUM", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        _print_output(bisect(model, w_x.jet[0], x_prev, f_tol, x_tol, mode=Solver.EXTREMUM))
                    elif n == 2:
                        _print_output(newton(model, w_x.jet[0], f_tol, x_tol, mode=Solver.EXTREMUM))
                    else:
                        _print_output(householder(model, w_x.jet[0], n, f_tol, x_tol, mode=Solver.EXTREMUM))
                # noinspection PyUnboundLocalVariable
                if f_dash_dash_prev * w_f.jet[2] < 0.0:
                    if f_dash_dash_prev > w_f.jet[2]:
                        print("+INFLECTION", file=stderr)
                    else:
                        print("-INFLECTION", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        _print_output(bisect(model, w_x.jet[0], x_prev, f_tol, x_tol, mode=Solver.INFLECTION))
                    elif n == 2:
                        _print_output(newton(model, w_x.jet[0], f_tol, x_tol, mode=Solver.INFLECTION))
                    else:
                        _print_output(householder(model, w_x.jet[0], n, f_tol, x_tol, mode=Solver.INFLECTION))
        x_prev = w_x.jet[0]
        f_prev = w_f.jet[0]
        f_dash_prev = w_f.jet[1]
        f_dash_dash_prev = w_f.jet[2]
