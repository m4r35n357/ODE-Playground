#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum
from ad import Series, to_mpfr


class Sense(Enum):
    POSITIVE = '+'
    NEGATIVE = '-'


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['count', 'sense', 'mode', 'x', 'f', 'dx'])


def bisect(model, ax, bx, f_tol, x_tol, max_it, sense, target, mode):
    a = Series.get(3, ax, variable=True)
    b = Series.get(3, bx, variable=True)
    c = Series.get(3)
    fc = Series.get(3, 1)
    f_sign = ~ model(a, target)
    delta = 1
    counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / 2
        fc = ~ model(c, target)
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        delta = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + target, dx=delta)


def newton(model, initial, f_tol, x_tol, max_it, sense, target, mode):
    x = Series.get(2 + mode.value, initial, variable=True)
    f = Series.get(2 + mode.value, 1)
    delta = 1
    counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        f = ~ model(x, target)
        delta = - f.jet[mode.value] / f.jet[1 + mode.value]
        x.val += delta
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + target, dx=delta)


def analyze(model, mode, x0, x1, steps, f_tol, x_tol, max_it, order):
    x_prev = f_prev = f_dash_prev = f_dash_dash_prev = result = None
    w_x = Series.get(order, x0, variable=True)
    target = to_mpfr(0.0)
    if mode != 0:
        if mode == 1:
            print("Bisection", file=stderr)
        elif mode == 2:
            print("Newton's method", file=stderr)
    for k in range(steps):
        w_x.val = x0 + k * (x1 - x0) / steps
        w_f = ~ model(w_x, target)
        print("{:.6e} {}".format(w_x.val, w_f))
        if mode != 0:
            if k > 0:
                if f_prev * w_f.val < 0.0:
                    if f_prev > w_f.val:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, target, Solver.ROOT)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, target, Solver.ROOT)
                    yield result
                if f_dash_prev * w_f.jet[1] < 0.0:
                    if f_dash_prev > w_f.jet[1]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, target, Solver.EXTREMUM)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, target, Solver.EXTREMUM)
                    yield result
                if f_dash_dash_prev * w_f.jet[2] < 0.0:
                    if f_dash_dash_prev > w_f.jet[2]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, target, Solver.INFLECTION)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, target, Solver.INFLECTION)
                    yield result
        else:
            yield w_x.val, w_f.jet
        x_prev = w_x.val
        f_prev = w_f.val
        f_dash_prev = w_f.jet[1]
        f_dash_dash_prev = w_f.jet[2]


print(__name__ + " module loaded", file=stderr)
