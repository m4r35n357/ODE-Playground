#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from ad import Series


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
    FLAT = '_'
    INCREASING = '/'
    DECREASING = '\\'


@unique
class Mode(Enum):
    ROOT = 0
    MIN_MAX = 1
    INFLECT = 2


class Bracketed(namedtuple('BracketedType', ['method', 'a', 'b', 'f', 'δx', 'count', 'sense', 'mode'])):
    def __str__(self):
        return f"{self.method} a={self.a:+.15e} b={self.b:+.15e} f={self.f:+.15e} δx={self.δx:+.15e} {self.sense}{self.mode} count={self.count}"


class Derivative(namedtuple('DerivativeType', ['method', 'x', 'f', 'δx', 'count', 'sense', 'mode'])):
    def __str__(self):
        return f"{self.method} x={self.x:+.15e} f={self.f:+.15e} δx={self.δx:+.15e} {self.sense}{self.mode} count={self.count}"


def bisect(model, xa, xb, εf=1e-15, εx=1e-15, limit=1001, sense=Sense.FLAT, y=0.0, mode=Mode.ROOT, debug=False):
    method = Solver.BI
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fc = Series.get(3, 1)
    f_sign = ~model(a) - y
    δx = counter = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a + b) / 2
        fc = ~model(c) - y
        if f_sign.jet[mode.value] * fc.jet[mode.value] < 0.0:
            b = c
        else:
            a = c
        δx = b.jet[mode.value] - a.jet[mode.value]
        if debug:
            print(Bracketed(method=method.name, count=counter, sense=sense.value, mode=mode.name, a=a.val, b=b.val, f=fc.val, δx=δx),
                  file=stderr)
        counter += 1
        if counter == limit:
            break
    return Bracketed(method=method.name, count=counter - 1, sense=sense.value, mode=mode.name, a=a.val, b=b.val, f=fc.val, δx=δx)


def falsi(model, xa, xb, εf=1e-15, εx=1e-15, limit=1001, sense=Sense.FLAT, y=0.0, mode=Mode.ROOT, ill=True, debug=False):
    method = Solver.FI if ill else Solver.FP
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~model(a) - y, ~model(b) - y, Series.get(3, 1)
    δx = counter = 1
    side = 0
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~model(c) - y
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
        δx = b.jet[mode.value] - a.jet[mode.value]
        if debug:
            print(Bracketed(method=method.name, count=counter, sense=sense.value, mode=mode.name, a=a.val, b=b.val, f=fc.val, δx=δx),
                  file=stderr)
        counter += 1
        if counter == limit:
            break
    return Bracketed(method=method.name, count=counter - 1, sense=sense.value, mode=mode.name, a=a.val, b=b.val, f=fc.val, δx=δx)


def secant(model, xa, xb, εf=1e-15, εx=1e-15, limit=1001, sense=Sense.FLAT, y=0.0, mode=Mode.ROOT, debug=False):
    method = Solver.SC
    a, b, c = Series.get(3, xa), Series.get(3, xb), Series.get(3)
    fa, fb, fc = ~model(a) - y, ~model(b) - y, Series.get(3, 1)
    δx = counter = 1
    while abs(fc.jet[mode.value]) > εf or abs(δx) > εx:
        c = (a * fb - b * fa) / (fb - fa)
        fc = ~model(c) - y
        b, fb = a, fa
        a, fa = c, fc
        δx = b.jet[mode.value] - a.jet[mode.value]
        if debug:
            print(Derivative(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=c.val, f=fc.val, δx=δx),
                  file=stderr)
        counter += 1
        if counter == limit:
            break
    return Derivative(method=method.name, count=counter - 1, sense=sense.value, mode=mode.name, x=c.val, f=fc.val, δx=δx)


def newton(model, x0, εf=1e-15, εx=1e-15, limit=1001, sense=Sense.FLAT, y=0.0, mode=Mode.ROOT, debug=False):
    method = Solver.NT
    x = Series.get(2 + mode.value, x0).var  # make x variable for AD
    f = Series.get(2 + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = ~model(x) - y
        δx = - f.jet[mode.value] / f.jet[1 + mode.value]
        x += δx
        if debug:
            print(Derivative(method=method.name, count=counter, sense=sense.value, mode=mode.name,x=x.val, f=f.val, δx=δx),
                  file=stderr)
        counter += 1
        if counter == limit:
            break
    return Derivative(method=method.name, count=counter - 1, sense=sense.value, mode=mode.name, x=x.val, f=f.val, δx=δx)


def householder(model, initial, n, εf=1e-15, εx=1e-15, limit=1001, sense=Sense.FLAT, y=0.0, mode=Mode.ROOT, debug=False):
    method = Solver.H1 if n == 2 else (Solver.H2 if n == 3 else (Solver.H3 if n == 4 else (Solver.H4 if n == 5 else Solver.NA)))
    x = Series.get(n + mode.value, initial).var  # make x variable for AD
    f = Series.get(n + mode.value, 1)
    δx = counter = 1
    while abs(f.jet[mode.value]) > εf or abs(δx) > εx:
        f = model(x) - y
        r = ~(1 / f)
        δx = r.jet[n - 2 + mode.value] / r.jet[n - 1 + mode.value]
        x += δx * (n - 1 + mode.value)
        if debug:
            print(Derivative(method=method.name, count=counter, sense=sense.value, mode=mode.name, x=x.val, f=f.val, δx=δx),
                  file=stderr)
        counter += 1
        if counter == limit:
            break
    return Derivative(method=method.name, count=counter - 1, sense=sense.value, mode=mode.name, x=x.val, f=f.val, δx=δx)


def analyze(model, method, x0, x1, steps, εf, εx, limit, order):
    x_prev = f0_prev = f1_prev = f2_prev = None
    for k in range(steps):
        step = (x1 - x0) / (steps - 1)
        x = Series.get(order, x0 + k * step).var  # make x a variable to see derivatives!
        f = ~model(x)
        print(f"{x.val:.6e} {f}")
        if method != Solver.NA:
            if k > 0:
                if f0_prev * f.val <= 0.0:
                    s = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, ill=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s)
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
                if f1_prev * f.jet[1] <= 0.0:
                    s = Sense.DECREASING if f1_prev > f.jet[1] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX, ill=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.MIN_MAX)
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
                if f2_prev * f.jet[2] <= 0.0:
                    s = Sense.DECREASING if f2_prev > f.jet[2] else Sense.INCREASING
                    if method == Solver.BI:
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.FP:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT, ill=False)
                    elif method == Solver.FI:
                        yield falsi(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
                    elif method == Solver.SC:
                        yield secant(model, x.val, x.val + step, εf=εf, εx=εx, limit=limit, sense=s, mode=Mode.INFLECT)
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


print(__name__ + " module loaded", file=stderr)
