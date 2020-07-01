#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from ad import Context, Series, Dual

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

def s_bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    a, b, c = Series.get(3, xa).var, Series.get(3, xb).var, Series.get(3)
    f_sign = model(a).jet[mode.value]
    fc = δx = i = 1
    while i <= limit and (abs(fc) > εf or abs(δx) > εx):
        c = 0.5 * (a + b)
        fc = (~ model(c)).jet[mode.value]
        if f_sign * fc < 0.0:
            b = c
        elif f_sign * fc > 0.0:
            a = c
        else:
            break
        δx = b.val - a.val
        i += 1
        if debug:
            print(Result(method=Solver.BI.name, count=i-1, sense=sense.value, mode=mode.name, x=c.val, f=fc, δx=δx), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, sense=sense.value, mode=mode.name, x=c.val, f=fc, δx=δx)

def s_newton(model, x0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, mode=Mode.ROOT___, debug=False):
    x, f = Series.get(2 + mode.value, x0).var, [1.0, 0.0]
    δx = i = 1
    while i <= limit and (abs(f[0]) > εf or abs(δx) > εx):
        f = (~ model(x)).jet[mode.value : 2 + mode.value]
        δx = - f[0] / f[1]
        x += δx
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, sense=sense.value, mode=mode.name, x=x.val, f=f[0], δx=δx), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, sense=sense.value, mode=mode.name, x=x.val, f=f[0], δx=δx)

def s_analyze(model, method, x0, x1, steps, εf, εx, limit, mode=Mode.ALL, console=True, debug=False):
    x_prev = f0_prev = f1_prev = f2_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        x = Series.get(3, x0 + k * step).var
        f = ~ model(x)
        if not console:
            print(f'{x.val:.{Context.places}e} {f}')
        if method != Solver.NA:
            if k > 0:
                if (mode == Mode.ROOT___ or mode == Mode.ALL) and f0_prev * f.val < 0.0:
                    sense = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Solver.BI:
                        yield s_bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=sense, debug=debug)
                    if method == Solver.NT:
                        yield s_newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense, debug=debug)
                if (mode == Mode.MIN_MAX or mode == Mode.ALL) and f1_prev * f.jet[Mode.MIN_MAX.value] < 0.0:
                    sense = Sense.DECREASING if f1_prev > f.jet[Mode.MIN_MAX.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield s_bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=sense, mode=Mode.MIN_MAX, debug=debug)
                    if method == Solver.NT:
                        yield s_newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense, mode=Mode.MIN_MAX, debug=debug)
                if (mode == Mode.INFLECT or mode == Mode.ALL) and f2_prev * f.jet[Mode.INFLECT.value] < 0.0:
                    sense = Sense.DECREASING if f2_prev > f.jet[Mode.INFLECT.value] else Sense.INCREASING
                    if method == Solver.BI:
                        yield s_bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=sense, mode=Mode.INFLECT, debug=debug)
                    if method == Solver.NT:
                        yield s_newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense, mode=Mode.INFLECT, debug=debug)
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[Mode.MIN_MAX.value]
        f2_prev = f.jet[Mode.INFLECT.value]

def d_bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    a, b, c = Dual.get(xa), Dual.get(xb), Dual.get(3)
    f_sign = model(Dual.get(xa)).val
    δx = fc = i = 1
    while i <= limit and abs(fc) > εf or abs(δx) > εx:
        c = 0.5 * (a + b)
        fc = model(c).val
        if f_sign * fc < 0.0:
            b = c
        elif f_sign * fc > 0.0:
            a = c
        else:
            break
        δx = b.val - a.val
        i += 1
        if debug:
            print(Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=fc, δx=δx, mode=Mode.ROOT___), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=fc, δx=δx, mode=Mode.ROOT___)

def d_newton(model, x0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    x, f = Dual.get(x0).var, Dual.get(1)
    δx = i = 1
    while i <= limit and abs(f.val) > εf or abs(δx) > εx:
        f = model(x)
        δx = - f.val / f.der
        x += δx
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=x.val, f=f.val, δx=δx, mode=Mode.ROOT___), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=x.val, f=f.val, δx=δx, mode=Mode.ROOT___)

def d_analyze(model, method, x0, x1, steps, εf, εx, limit, console=True):
    x_prev = f0_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        x = Dual.get(x0 + k * step).var
        f = model(x)
        if not console:
            print(f'{x.val:.{Context.places}e} {f}')
        if method != Solver.NA:
            if k > 0:
                if f0_prev * f.val < 0.0:
                    sense = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Solver.BI:
                        yield d_bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=sense)
                    if method == Solver.NT:
                        yield d_newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense)
        x_prev = x.val
        f0_prev = f.val

print(f'{__name__} module loaded', file=stderr)
