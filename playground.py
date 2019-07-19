#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum
from gmpy2 import zero
from ad import Series, to_mpfr
from functions import x_step


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


def analyze(model, mode, x0, x1, steps, f_tol, x_tol, max_it, order):
    x_prev = f0_prev = f1_prev = f2_prev = result = None
    if mode != 0:
        if mode == 1:
            print("Bisection", file=stderr)
        elif mode == 2:
            print("Newton's method", file=stderr)
    for k in range(steps):
        w_x = Series.get(order, x_step(x0, x1, steps, k)).var  # make x a variable to see derivatives!
        w_f = ~model(w_x)
        print(f"{w_x.val:.6e} {w_f}")
        if mode != 0:
            if k > 0:
                if f0_prev * w_f.val <= zero(+1):
                    if f0_prev > w_f.val:
                        sense = Sense.DECREASING
                    else:
                        sense = Sense.INCREASING
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense)
                    yield result
                if f1_prev * w_f.jet[1] <= zero(+1):
                    if f1_prev > w_f.jet[1]:
                        sense = Sense.DECREASING
                    else:
                        sense = Sense.INCREASING
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.EXTREMUM)
                    yield result
                if f2_prev * w_f.jet[2] <= zero(+1):
                    if f2_prev > w_f.jet[2]:
                        sense = Sense.DECREASING
                    else:
                        sense = Sense.INCREASING
                    if mode == 1:
                        result = bisect(model, w_x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif mode == 2:
                        result = newton(model, w_x.val, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    yield result
        else:
            yield w_x.val, w_f.jet
        x_prev = w_x.val
        f0_prev = w_f.val
        f1_prev = w_f.jet[1]
        f2_prev = w_f.jet[2]


print(__name__ + " module loaded", file=stderr)
