#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum
from taylor import t_jet, to_mpfr
from series import Series

# noinspection PyArgumentList
De_12 = to_mpfr(1.0e-12)


class Sense(Enum):
    POSITIVE = '+'
    NEGATIVE = '-'


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['count', 'sense', 'mode', 'x', 'f', 'dx'])


def bisect(model, ax, bx, f_tol, x_tol, max_it, sense, target=0, mode=Solver.ROOT):
    a = Series(t_jet(3, ax), variable=True)
    b = Series(t_jet(3, bx), variable=True)
    c = Series(t_jet(3))
    fc = Series(t_jet(3, 1))
    f_sign = ~ model(a, target)
    delta = 1
    counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / 2
        fc = ~ model(c, to_mpfr(target))
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        delta = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + to_mpfr(target), dx=delta)


def newton(model, initial, f_tol, x_tol, max_it, sense, target=0, mode=Solver.ROOT):
    x = Series(t_jet(2 + mode.value, initial), variable=True)
    f = Series(t_jet(2 + mode.value, 1))
    delta = 1
    counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        f = ~ model(x, to_mpfr(target))
        delta = - f.jet[mode.value] / f.jet[1 + mode.value]
        x.val += delta
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + to_mpfr(target), dx=delta)


def analyze(model, n, x0, x1, steps=1000, target=0, f_tol=De_12, x_tol=De_12, max_it=100, n_max=13):
    x_prev = f_prev = f_dash_prev = f_dash_dash_prev = result = None
    w_x = Series(t_jet(n_max, x0), variable=True)
    if n != 0:
        if n == 1:
            print("Bisection", file=stderr)
        elif n == 2:
            print("Newton's method", file=stderr)
    for k in range(steps):
        w_x.val = x0 + k * (x1 - x0) / steps
        w_f = ~ model(w_x, to_mpfr(target))
        print("{:.6e} {}".format(w_x.val, w_f))
        if n != 0:
            if k > 0:
                if f_prev * w_f.val < 0.0:
                    if f_prev > w_f.val:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, target=to_mpfr(target))
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, target=to_mpfr(target))
                    yield result
                if f_dash_prev * w_f.jet[1] < 0.0:
                    if f_dash_prev > w_f.jet[1]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    yield result
                if f_dash_dash_prev * w_f.jet[2] < 0.0:
                    if f_dash_dash_prev > w_f.jet[2]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    yield result
        else:
            yield w_x.val, w_f.jet
        x_prev = w_x.val
        f_prev = w_f.val
        f_dash_prev = w_f.jet[1]
        f_dash_dash_prev = w_f.jet[2]


print(__name__ + " module loaded", file=stderr)
