#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv, stderr
from collections import namedtuple
from enum import Enum
from gmpy2 import mpfr
from taylor import t_jet, D0, D1, D2
from series import Series


class Sense(Enum):
    POSITIVE = '+'
    NEGATIVE = '-'


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['count', 'sense', 'mode', 'x', 'f', 'dx'])


def bisect(model, ax, bx, f_tol, x_tol, max_it, sense, target=D0, mode=Solver.ROOT):
    a = Series(t_jet(3, ax), variable=True)
    b = Series(t_jet(3, bx), variable=True)
    c = Series(t_jet(3))
    fc = Series(t_jet(3, D1))
    f_sign = ~ model(a, target)
    delta = D1
    counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / D2
        fc = ~ model(c, target)
        if f_sign.jet[mode.value] * fc.jet[mode.value] < D0:
            b = c
        else:
            a = c
        delta = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + target, dx=delta)


def newton(model, initial, f_tol, x_tol, max_it, sense, target=D0, mode=Solver.ROOT):
    x = Series(t_jet(2 + mode.value, initial), variable=True)
    f = Series(t_jet(2 + mode.value, D1))
    delta = D1
    counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(delta) > x_tol:
        f = ~ model(x, target)
        delta = - f.jet[mode.value] / f.jet[1 + mode.value]
        x.val += delta
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + target, dx=delta)


def analyze(model, max_it):
    x_prev = f_prev = f_dash_prev = f_dash_dash_prev = result = None
    n_max = 13
    n = int(argv[1])
    assert n == 0 or n == 1 or n == 2
    # noinspection PyArgumentList
    x0 = mpfr(argv[2])
    # noinspection PyArgumentList
    x1 = mpfr(argv[3])
    assert x1 > x0
    steps = int(argv[4])
    assert steps > 0
    # noinspection PyArgumentList
    target = mpfr(argv[5])
    # noinspection PyArgumentList
    f_tol = mpfr(argv[6])
    # noinspection PyArgumentList
    x_tol = mpfr(argv[7])
    w_x = Series(t_jet(n_max, x0), variable=True)
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
                if f_prev * w_f.val < D0:
                    if f_prev > w_f.val:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, target=target)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, target=target)
                    yield result
                if f_dash_prev * w_f.jet[1] < D0:
                    if f_dash_prev > w_f.jet[1]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    yield result
                if f_dash_dash_prev * w_f.jet[2] < D0:
                    if f_dash_dash_prev > w_f.jet[2]:
                        sense = Sense.NEGATIVE
                    else:
                        sense = Sense.POSITIVE
                    if n == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif n == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    yield result
        x_prev = w_x.val
        f_prev = w_f.val
        f_dash_prev = w_f.jet[1]
        f_dash_dash_prev = w_f.jet[2]


print(__name__ + " module loaded", file=stderr)
