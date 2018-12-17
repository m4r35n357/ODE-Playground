#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv, stderr
from collections import namedtuple
from enum import Enum
from taylor import t_jet
from series import Series


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['count', 'mode', 'x', 'f', 'dx'])


def _print_output(output):
    print("    {:3d} {} {:22.15e} {:10.3e} {:10.3e}".format(*output), file=stderr)


def bisect(model, ax, bx, f_tol, x_tol, max_it, target=0.0, mode=Solver.ROOT):
    a = Series(t_jet(3, ax), diff=True)
    b = Series(t_jet(3, bx), diff=True)
    c = Series(t_jet(3))
    fc = Series(t_jet(3, 1.0))
    f_sign = ~ model(a, target)
    delta = 1.0
    counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / 2.0
        fc = ~ model(c, target)
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        delta = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter > max_it:
            break
    return Result(count=counter, mode=mode.name, x=c.val, f=fc.val + target, dx=delta)


def newton(model, initial, f_tol, x_tol, max_it, target=0.0, mode=Solver.ROOT):
    x = Series(t_jet(2 + mode.value, initial), diff=True)
    f = Series(t_jet(2 + mode.value, 1.0))
    delta = 1.0
    counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        f = ~ model(x, target)
        delta = - f.jet[mode.value] / f.jet[1 + mode.value]
        x.val += delta
        counter += 1
        if counter > max_it:
            break
    return Result(count=counter, mode=mode.name, x=x.val, f=f.val + target, dx=delta)


def analyze(model):
    n_max = 7
    max_it = 100
    n = int(argv[1])
    assert n == 0 or n == 1 or n == 2
    x0 = float(argv[2])
    x1 = float(argv[3])
    assert x1 > x0
    steps = int(argv[4])
    assert steps > 0
    target = float(argv[5])
    f_tol = float(argv[6])
    x_tol = float(argv[7])
    w_x = Series(t_jet(n_max, x0), diff=True)
    results = []
    if n != 0:
        if n == 1:
            print("Bisection", file=stderr)
        elif n == 2:
            print("Newton's method", file=stderr)
    for k in range(steps):
        w_x.val = x0 + k * (x1 - x0) / steps
        w_f = ~ model(w_x, target)
        print("{:.6e} {}".format(w_x.val, w_f))
        if n != 0:
            if k > 0:
                # noinspection PyUnboundLocalVariable
                if f_prev * w_f.val < 0.0:
                    print("ROOT", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, target=target)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, target=target)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
                # noinspection PyUnboundLocalVariable
                if f_dash_prev * w_f.jet[1] < 0.0:
                    if f_dash_prev > w_f.jet[1]:
                        print("MAXIMUM", file=stderr)
                    else:
                        print("MINIMUM", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, mode=Solver.EXTREMUM)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, mode=Solver.EXTREMUM)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
                # noinspection PyUnboundLocalVariable
                if f_dash_dash_prev * w_f.jet[2] < 0.0:
                    if f_dash_dash_prev > w_f.jet[2]:
                        print("+INFLECTION", file=stderr)
                    else:
                        print("-INFLECTION", file=stderr)
                    if n == 1:
                        # noinspection PyUnboundLocalVariable
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, mode=Solver.INFLECTION)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, mode=Solver.INFLECTION)
                        _print_output(result)
                        if result.count < max_it:
                            results.append(result)
        x_prev = w_x.val
        f_prev = w_f.val
        f_dash_prev = w_f.jet[1]
        f_dash_dash_prev = w_f.jet[2]
    return results


print(__name__ + " module loaded", file=stderr)
