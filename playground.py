#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from gmpy2 import zero
from ad import Series, to_mpfr


@unique
class Analysis(Enum):
    NA = "No analysis"
    BI = "Bisection method"
    FP = "False Position method"
    SC = "Secant method"
    NT = "Newton's method"
    H2 = "Householder's method (degree 2)"
    H3 = "Householder's method (degree 3)"


@unique
class Sense(Enum):
    FLAT = ''
    INCREASING = '+'
    DECREASING = '-'


@unique
class Solver(Enum):
    ROOT = 0
    MIN_MAX = 1
    INFLECTION = 2


Result = namedtuple('ResultType', ['method', 'count', 'sense', 'mode', 'x', 'f', 'δx'])


def bisect(model, xa, xb, f_tol=to_mpfr(1.0e-12), x_tol=to_mpfr(1.0e-12), max_it=1001, sense=Sense.FLAT, target=zero(+1), mode=Solver.ROOT):
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
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
    return Result(method=Analysis.BI.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + target, δx=δx)


def secant(model, xa, xb, f_tol=to_mpfr(1.0e-12), x_tol=to_mpfr(1.0e-12), max_it=1001, sense=Sense.FLAT, target=zero(+1), mode=Solver.ROOT, fp=False):
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~model(a) - target, ~model(b) - target, Series.get(3, 1)
    δx = counter = 1
    while abs(fc.jet[mode.value]) > f_tol or abs(δx) > x_tol:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~model(c) - target
        if abs(fc.jet[mode.value]) == zero(+1):
            break
        if fp:
            if fa.val * fc.val > zero(+1):
                a, fa = c, fc
            elif fb.val * fc.val > zero(+1):
                b, fb = c, fc
        else:
            b, fb = a, fa
            a, fa = c, fc
        δx = b.jet[mode.value] - a.jet[mode.value]
        counter += 1
        if counter == max_it:
            break
    method = Analysis.FP if fp else Analysis.SC
    return Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + target, δx=δx)


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
    return Result(method=Analysis.NT.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + target, δx=δx)


def householder(model, initial, n, f_tol=to_mpfr(1.0e-12), x_tol=to_mpfr(1.0e-12), max_it=100, sense=Sense.FLAT, target=zero(+1), mode=Solver.ROOT):
    x = Series.get(n + mode.value, initial).var
    f = Series.get(n + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > f_tol or abs(δx) > x_tol:
        f = model(x) - target
        r = ~ (1 / f)
        δx = r.jet[n - 2 + mode.value] / r.jet[n - 1 + mode.value]
        x += δx * (n - 1 + mode.value)
        counter += 1
        if counter == max_it:
            break
    method = Analysis.H2 if n == 3 else (Analysis.H3 if n == 4 else Analysis.NA)
    return Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + target, δx=δx)


def analyze(model, method, x0, x1, steps, f_tol, x_tol, max_it, order):
    x_prev = f0_prev = f1_prev = f2_prev = None
    print(method.value, file=stderr)
    for k in range(steps):
        step = (x1 - x0) / (steps - 1)
        x = Series.get(order, x0 + k * step).var  # make x a variable to see derivatives!
        f = ~model(x)
        print(f"{x.val:.6e} {f}")
        if method != Analysis.NA:
            if k > 0:
                if f0_prev * f.val <= zero(+1):
                    sense = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, f_tol, x_tol, max_it, sense)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, f_tol, x_tol, max_it, sense)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, f_tol, x_tol, max_it, sense)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, f_tol, x_tol, max_it, sense)
                if f1_prev * f.jet[1] <= zero(+1):
                    sense = Sense.DECREASING if f1_prev > f.jet[1] else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, f_tol, x_tol, max_it, sense, mode=Solver.MIN_MAX)
                if f2_prev * f.jet[2] <= zero(+1):
                    sense = Sense.DECREASING if f2_prev > f.jet[2] else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, f_tol, x_tol, max_it, sense, mode=Solver.INFLECTION)
        else:
            yield x.val, f.jet
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[1]
        f2_prev = f.jet[2]


print(__name__ + " module loaded", file=stderr)
