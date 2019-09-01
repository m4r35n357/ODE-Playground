#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from math import copysign
from ad import Context, Series

@unique
class Solver(Enum):
    NA = "No analysis, or all"
    BI = "Bisection method"
    NT = "Newton-Raphson method"

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
    ALL = 3

class Result(namedtuple('ResultType', ['method', 'x', 'f', 'δx', 'count', 'sense', 'mode'])):
    def __str__(self):
        return f'{self.method}  x: {self.x:+.{Context.places}e}  δx: {self.δx:+.{Context.places}e}  ' \
               f'f: {self.f:+.{Context.places}e}  {self.sense} {self.mode} {self.count}'

def bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.BI
    a, b, c = Series.get(3, xa).var, Series.get(3, xb).var, Series.get(3)
    fc = Series.get(3, 1)
    f_sign = copysign(1, (~ model(a)).jet[mode.value])
    δx = count = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a + b) / 2
        fc = ~ model(c)
        if f_sign * fc.jet[mode.value] < 0.0:
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

def newton(model, x0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    m = Solver.NT
    x = Series.get(2 + mode.value, x0).var  # make x variable for AD
    f = Series.get(2 + mode.value, 1)
    δx = count = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = ~ model(x)
        δx = - f.jet[mode.value] / f.jet[1 + mode.value]
        x += δx
        if debug:
            print(Result(method=m.name, count=count, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx))
        count += 1
        if count == limit + 1:
            break
    return Result(method=m.name, count=count - 1, sense=sense.value, mode=mode.name, x=x.val, f=f.jet[mode.value], δx=δx)

def analyze(model, method, x0, x1, steps, εf, εx, limit, order, mode=Mode.ALL, console=True):
    x_prev = f0_prev = f1_prev = f2_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        x = Series.get(order, x0 + k * step).var  # make x a variable to see derivatives!
        f = ~ model(x)
        if not console:
            print(f'{x.val:.{Context.places}e} {f}')
        if method != Solver.NA:
            if k > 0:
                if (mode == Mode.ROOT___ or mode == Mode.ALL) and f0_prev * f.val < 0.0:
                    s = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s)
                    if method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s)
                if (mode == Mode.MIN_MAX or mode == Mode.ALL) and f1_prev * f.jet[Mode.MIN_MAX.value] < 0.0:
                    s = Sense.DECREASING if f1_prev > f.jet[Mode.MIN_MAX.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    if method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                if (mode == Mode.INFLECT or mode == Mode.ALL) and f2_prev * f.jet[Mode.INFLECT.value] < 0.0:
                    s = Sense.DECREASING if f2_prev > f.jet[Mode.INFLECT.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    if method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[Mode.MIN_MAX.value]
        f2_prev = f.jet[Mode.INFLECT.value]

print(f'{__name__} module loaded', file=stderr)
