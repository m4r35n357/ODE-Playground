#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from ad import Series


@unique
class Analysis(Enum):
    NA = "No analysis, or all"
    BI = "Bisection method"
    FP = "False Position method (Illinois Algorithm)"
    SC = "Secant method"
    NT = "Newton-Raphson method"
    H1 = "Householder's method, degree 1 (Newton)"
    H2 = "Householder's method, degree 2 (Halley)"
    H3 = "Householder's method, degree 3"
    H4 = "Householder's method, degree 4"


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


Result = namedtuple('ResultType', ['method', 'x', 'f', 'δx', 'count', 'sense', 'mode'])


def bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=1001, sense=Sense.FLAT, y=0.0, mode=Solver.ROOT, debug=False):
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fc = Series.get(3, 1)
    f_sign = ~model(a) - y
    δx = counter = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a + b) / 2
        fc = ~model(c) - y
        if abs(fc.jet[mode.value]) == 0.0:
            break
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        δx = b.jet[mode.value] - a.jet[mode.value]
        if debug:
            print(Result(method=Analysis.BI.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + y, δx=δx))
        counter += 1
        if counter == limit:
            break
    return Result(method=Analysis.BI.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + y, δx=δx)


def secant(model, xa, xb, εf=1e-12, εx=1e-12, limit=1001, sense=Sense.FLAT, y=0.0, mode=Solver.ROOT, fp=False, ill=True, debug=False):
    method = Analysis.FP if fp else Analysis.SC
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~model(a) - y, ~model(b) - y, Series.get(3, 1)
    δx = counter = 1
    side = 0
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~model(c) - y
        if abs(fc.jet[mode.value]) == 0.0:
            break
        if fp:
            if fa.val * fc.val > 0.0:
                a, fa = c, fc
                if ill:
                    if side == -1:
                        fb *= 0.5
                    side = -1
            elif fb.val * fc.val > 0.0:
                b, fb = c, fc
                if ill:
                    if side == +1:
                        fa *= 0.5
                    side = +1
        else:
            b, fb = a, fa
            a, fa = c, fc
        δx = b.jet[mode.value] - a.jet[mode.value]
        if debug:
            print(Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + y, δx=δx))
        counter += 1
        if counter == limit:
            break
    return Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val + y, δx=δx)


def newton(model, x0, εf=1e-12, εx=1e-12, limit=1001, sense=Sense.FLAT, y=0.0, mode=Solver.ROOT, debug=False):
    x = Series.get(2 + mode.value, x0).var  # make x a variable to use derivatives!
    f = Series.get(2 + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = ~model(x) - y
        if abs(f.jet[mode.value]) == 0.0:
            break
        δx = - f.jet[mode.value] / f.jet[1 + mode.value]
        x += δx
        if debug:
            print(Result(method=Analysis.NT.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + y, δx=δx))
        counter += 1
        if counter == limit:
            break
    return Result(method=Analysis.NT.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + y, δx=δx)


def householder(model, initial, n, εf=1e-12, εx=1e-12, limit=100, sense=Sense.FLAT, y=0.0, mode=Solver.ROOT, debug=False):
    method = Analysis.H1 if n == 2 else (Analysis.H2 if n == 3 else (Analysis.H3 if n == 4 else (Analysis.H4 if n == 5 else Analysis.NA)))
    x = Series.get(n + mode.value, initial).var
    f = Series.get(n + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = model(x) - y
        if abs(f.jet[mode.value]) == 0.0:
            break
        r = ~(1 / f)
        δx = r.jet[n - 2 + mode.value] / r.jet[n - 1 + mode.value]
        x += δx * (n - 1 + mode.value)
        if debug:
            print(Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + y, δx=δx))
        counter += 1
        if counter == limit:
            break
    return Result(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val + y, δx=δx)


def analyze(model, method, x0, x1, steps, εf, εx, limit, order):
    x_prev = f0_prev = f1_prev = f2_prev = None
    print(method.value, file=stderr)
    for k in range(steps):
        step = (x1 - x0) / (steps - 1)
        x = Series.get(order, x0 + k * step).var  # make x a variable to see derivatives!
        f = ~model(x)
        print(f"{x.val:.6e} {f}")
        if method != Analysis.NA:
            if k > 0:
                if f0_prev * f.val <= 0.0:
                    sense = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, εf, εx, limit, sense)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, εf, εx, limit, sense)
                    elif method == Analysis.H1:
                        yield householder(model, x.val, 2, εf, εx, limit, sense)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, εf, εx, limit, sense)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, εf, εx, limit, sense)
                    elif method == Analysis.H4:
                        yield householder(model, x.val, 5, εf, εx, limit, sense)
                if f1_prev * f.jet[1] <= 0.0:
                    sense = Sense.DECREASING if f1_prev > f.jet[1] else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense, mode=Solver.MIN_MAX, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H1:
                        yield householder(model, x.val, 2, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                    elif method == Analysis.H4:
                        yield householder(model, x.val, 5, εf, εx, limit, sense, mode=Solver.MIN_MAX)
                if f2_prev * f.jet[2] <= 0.0:
                    sense = Sense.DECREASING if f2_prev > f.jet[2] else Sense.INCREASING
                    if method == Analysis.BI:
                        yield bisect(model, x.val, x_prev, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.FP:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense, mode=Solver.INFLECTION, fp=True)
                    elif method == Analysis.SC:
                        yield secant(model, x.val, x.val + step, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.NT:
                        yield newton(model, x.val, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H1:
                        yield householder(model, x.val, 2, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H2:
                        yield householder(model, x.val, 3, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H3:
                        yield householder(model, x.val, 4, εf, εx, limit, sense, mode=Solver.INFLECTION)
                    elif method == Analysis.H4:
                        yield householder(model, x.val, 5, εf, εx, limit, sense, mode=Solver.INFLECTION)
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[1]
        f2_prev = f.jet[2]


print(__name__ + " module loaded", file=stderr)
