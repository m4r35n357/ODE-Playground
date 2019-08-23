#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from ad import Context, Series


@unique
class Solver(Enum):
    NA = "No analysis, or all"
    BI = "Bisection method"
    FP = "False Position method"
    FI = "False Position method (Illinois Algorithm)"
    SC = "Secant method"
    NT = "Newton-Raphson method"
    H1 = "Householder's method, degree 1 (Newton)"
    H2 = "Householder's method, degree 2 (Halley)"
    H3 = "Householder's method, degree 3"
    H4 = "Householder's method, degree 4"


@unique
class Sense(Enum):
    INCREASING = '/'
    DECREASING = '\\'
    FLAT = '_'


@unique
class Mode(Enum):
    ROOT___ = 0
    MIN_MAX = 1
    INFLECT = 2


class Result(namedtuple('ResultType', ['method', 'x', 'f', 'δx', 'count', 'sense', 'mode'])):
    def __str__(self):
        return f'{self.method}  x: {self.x:+.{Context.places}e}  δx: {self.δx:+.{Context.places}e}  ' \
               f'f: {self.f:+.{Context.places}e}  {self.sense} {self.mode} {self.count}'


def bisect(model, xa, xb, y=0.0, εf=1e-15, εx=1e-15, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.BI
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fc = Series.get(3, 1)
    f_sign = ~ model(a) - y
    δx = count = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a + b) / 2
        fc = ~ model(c) - y
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        δx = b.val - a.val
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx)


def falsi(model, xa, xb, y=0.0, εf=1e-15, εx=1e-15, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, illinois=True, debug=False):
    m = Solver.FI if illinois else Solver.FP
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~ model(a) - y, ~ model(b) - y, Series.get(3, 1)
    δx = count = 1
    side = 0
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~ model(c) - y
        if fa.jet[mode.value] * fc.jet[mode.value] > 0.0:
            a, fa = c, fc
            if illinois:
                if side == -1:
                    fb *= 0.5
                side = -1
        elif fb.jet[mode.value] * fc.jet[mode.value] > 0.0:
            b, fb = c, fc
            if illinois:
                if side == +1:
                    fa *= 0.5
                side = +1
        δx = b.val - a.val
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx)


def secant(model, xa, xb, y=0.0, εf=1e-15, εx=1e-15, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.SC
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~ model(a) - y, ~ model(b) - y, Series.get(3, 1)
    δx = count = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~ model(c) - y
        b, fb = a, fa
        a, fa = c, fc
        δx = b.val - a.val
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=c.val, f=fc.jet[mode.value], δx=δx)


def newton(model, x0, y=0.0, εf=1e-15, εx=1e-15, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.NT
    x = Series.get(2 + mode.value, x0).var  # make x variable for AD
    f = Series.get(2 + mode.value, 1)
    δx = count = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = ~ model(x) - y
        δx = - f.jet[mode.value] / f.jet[1 + mode.value]
        x += δx
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx)


def householder(model, x0, n, y=0.0, εf=1e-15, εx=1e-15, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.H1 if n == 2 else (Solver.H2 if n == 3 else (Solver.H3 if n == 4 else (Solver.H4 if n == 5 else Solver.NA)))
    x = Series.get(n + mode.value, x0).var  # make x variable for AD
    f = Series.get(n + mode.value, 1)
    δx = count = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = model(x) - y
        r = ~ (1 / f)
        δx = r.jet[n - 2 + mode.value] / r.jet[n - 1 + mode.value]
        x += δx * (n - 1 + mode.value)
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx)


def analyze(model, method, x0, x1, steps, εf, εx, limit, order):
    x_prev = f0_prev = f1_prev = f2_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        x = Series.get(order, x0 + k * step).var  # make x a variable to see derivatives!
        f = ~ model(x)
        print(f'{x.val:.{Context.places}e} {f}')
        if method != Solver.NA:
            if k > 0:
                if f0_prev * f.val < 0.0:
                    s = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, illinois=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.H1:
                        yield householder(model, x.val, 2, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.H2:
                        yield householder(model, x.val, 3, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.H3:
                        yield householder(model, x.val, 4, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.H4:
                        yield householder(model, x.val, 5, εf=εf, εx=εx, limit=limit, sense=s)
                if f1_prev * f.jet[Mode.MIN_MAX.value] < 0.0:
                    s = Sense.DECREASING if f1_prev > f.jet[Mode.MIN_MAX.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX, illinois=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.H1:
                        yield householder(model, x.val, 2, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.H2:
                        yield householder(model, x.val, 3, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.H3:
                        yield householder(model, x.val, 4, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.H4:
                        yield householder(model, x.val, 5, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                if f2_prev * f.jet[Mode.INFLECT.value] < 0.0:
                    s = Sense.DECREASING if f2_prev > f.jet[Mode.INFLECT.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT, illinois=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.H1:
                        yield householder(model, x.val, 2, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.H2:
                        yield householder(model, x.val, 3, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.H3:
                        yield householder(model, x.val, 4, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.H4:
                        yield householder(model, x.val, 5, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[1]
        f2_prev = f.jet[2]


print(f'{__name__} module loaded', file=stderr)
