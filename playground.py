#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum
from gmpy2 import zero
from ad import Series, to_mpfr
from functions import x_step


class Analysis(Enum):
    SKIP = 0
    BISECTION = 1
    NEWTON = 2


class Sense(Enum):
    FLAT = ''
    INCREASING = '+'
    DECREASING = '-'


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['count', 'sense', 'mode', 'x', 'f', 'δx'])


def bisect(model, xa, xb, f_tol=to_mpfr(1.0e-12), x_tol=to_mpfr(1.0e-12), max_it=1001, sense=Sense.FLAT, target=zero(+1), mode=Solver.ROOT):
    a = Series.get(3, xa)
    b = Series.get(3, xb)
    c = Series.get(3)
    fc = Series.get(3, 1)
    f_sign = ~model(a) - target
    δx = counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(δx) > x_tol:
        c = (a + b) / 2
        fc = ~model(c) - target
        if abs(fc.jet[mode.value]) == zero(+1):
            break
        if f_sign.jet[mode.value] * fc.jet[mode.value] < zero(+1):
            b = c
        else:
            a = c
        δx = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + target, δx=δx)


def newton(model, x0, f_tol=to_mpfr(1.0e-12), x_tol=to_mpfr(1.0e-12), max_it=1001, sense=Sense.FLAT, target=zero(+1), mode=Solver.ROOT):
    x = Series.get(2 + mode.value, x0).var  # make x a variable to use derivatives!
    f = Series.get(2 + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(δx) > x_tol:
        f = ~model(x) - target
        if abs(f.jet[mode.value]) == zero(+1):
            break
        δx = - f.jet[mode.value] / f.jet[1 + mode.value]
        x += δx
        counter += 1
        if counter == max_it:
            break
    return Result(count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + target, δx=δx)


def analyze(model, method, x0, x1, steps, f_tol, x_tol, max_it, order):
    x_prev = f0_prev = f1_prev = f2_prev = None
    if method != Analysis.SKIP.value:
        if method == Analysis.BISECTION.value:
            print("Bisection method", file=stderr)
        elif method == Analysis.NEWTON.value:
            print("Newton's method", file=stderr)
    for k in range(steps):
        w_x = Series.get(order, x_step(x0, x1, steps, k)).var  # make x a variable to see derivatives!
        w_f = ~model(w_x)
        print(f"{w_x.val:.6e} {w_f}")
        if method != Analysis.SKIP.value:
            if k > 0:
                if f0_prev * w_f.val <= zero(+1):
                    sense = Sense.DECREASING if f0_prev > w_f.val else Sense.INCREASING
                    if method == Analysis.BISECTION.value:
                        yield bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense)
                    elif method == Analysis.NEWTON.value:
                        yield newton(model, w_x.val, f_tol, x_tol, max_it, sense)
                if f1_prev * w_f.jet[1] <= zero(+1):
                    sense = Sense.DECREASING if f1_prev > w_f.jet[1] else Sense.INCREASING
                    if method == Analysis.BISECTION.value:
                        yield bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    elif method == Analysis.NEWTON.value:
                        yield newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                if f2_prev * w_f.jet[2] <= zero(+1):
                    sense = Sense.DECREASING if f2_prev > w_f.jet[2] else Sense.INCREASING
                    if method == Analysis.BISECTION.value:
                        yield bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.NEWTON.value:
                        yield newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
        else:
            yield w_x.val, w_f.jet
        x_prev = w_x.val
        f0_prev = w_f.val
        f1_prev = w_f.jet[1]
        f2_prev = w_f.jet[2]


print(__name__ + " module loaded", file=stderr)
