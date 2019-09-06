#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from math import sin, cos, sinh, cosh, tan, tanh, exp, log, copysign
from matplotlib import pyplot

class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value=0.0):
        return cls(value if isinstance(value, float) else float(value), 0.0)

    def __str__(self):
        return f'{self.val:+.{Context.places}e} {self.der:+.{Context.places}e}'

    def __abs__(self):
        return Dual(abs(self.val), self.der if self.val > 0.0 else (- self.der if self.val < 0.0 else 0.0))

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val + o.val, self.der + o.der)
        elif isinstance(o, (float, int)):
            return Dual(self.val + o, self.der)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __radd__(self, o):
        return self.__add__(o)

    def __sub__(self, o):
        return self.__add__(- o)

    def __rsub__(self, o):
        return (- self).__add__(o)

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val * o.val, self.der * o.val + self.val * o.der)
        elif isinstance(o, (float, int)):
            return Dual(self.val * o, self.der * o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self.__mul__(o)

    def __truediv__(self, o):
        if isinstance(o, Dual):
            assert abs(o.val) != 0.0, f"other.val = {o.val}"
            return Dual(self.val / o.val, (self.der * o.val - self.val * o.der) / o.val**2)
        elif isinstance(o, (float, int)):
            assert abs(o) != 0.0, f"other = {o}"
            return self.__mul__(1.0 / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        return Dual(o, 0.0).__truediv__(self)

    def __pow__(self, o):
        if isinstance(o, int):
            i_pow = Dual(self.val, self.der)
            for _ in range(abs(o) - 1):
                i_pow = i_pow.__mul__(self)
            return i_pow if o > 0 else (Dual(1.0, 0.0).__truediv__(i_pow) if o < 0 else Dual(1.0, 0.0))
        else:
            if isinstance(o, Dual):
                return (self.ln.__mul__(o)).exp
            elif isinstance(o, float):
                assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
                return Dual(self.val**o, o * self.val**(o - 1) * self.der)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self.__mul__(log(o))).exp

    @property
    def exp(self):
        exp_val = exp(self.val)
        return Dual(exp_val, self.der * exp_val)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return Dual(log(self.val), self.der / self.val)

    @property
    def sin(self):
        return Dual(sin(self.val), self.der * cos(self.val))

    @property
    def cos(self):
        return Dual(cos(self.val), - self.der * sin(self.val))

    @property
    def tan(self):
        return Dual(tan(self.val), self.der * (1.0 + tan(self.val)**2))

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.der * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.der * sinh(self.val))

    @property
    def tanh(self):
        return Dual(tanh(self.val), self.der * (1.0 - tanh(self.val)**2))

    @property
    def var(self):
        return Dual(self.val, 1.0)

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

class Result(namedtuple('ResultType', ['method', 'x', 'f', 'δx', 'count', 'sense'])):
    def __str__(self):
        return f'{self.method}  x: {self.x:+.{Context.places}e}  δx: {self.δx:+.{Context.places}e}  ' \
               f'f: {self.f:+.{Context.places}e}  {self.sense} {self.count}'

def bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    a, b, c = Dual.get(xa).var, Dual.get(xb).var, Dual.get(3)
    fc = Dual.get(1)
    f_sign = copysign(1, model(a).val)
    δx = count = 1
    while abs(fc.val) > εf or abs(δx) > εx:
        c = (a + b) / 2
        fc = model(c)
        if f_sign * fc.val < 0.0:
            b = c
        elif f_sign * fc.val > 0.0:
            a = c
        else:
            break
        δx = b.val - a.val
        if debug:
            print(Result(method=Solver.BI.name, count=count, sense=sense.value, x=c.val, f=fc.val, δx=δx), file=stderr)
        count += 1
        if count == limit + 1:
            break
    return Result(method=Solver.BI.name, count=count - 1, sense=sense.value, x=c.val, f=fc.val, δx=δx)

def newton(model, x0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    x = Dual.get(x0).var
    f = Dual.get(1)
    δx = count = 1
    while abs(f.val) > εf or abs(δx) > εx:
        f = model(x)
        δx = - f.val / f.der
        x += δx
        if debug:
            print(Result(method=Solver.NT.name, count=count, sense=sense.value, x=x.val, f=f.val, δx=δx), file=stderr)
        count += 1
        if count == limit + 1:
            break
    return Result(method=Solver.NT.name, count=count - 1, sense=sense.value, x=x.val, f=f.val, δx=δx)

def analyze(model, method, x0, x1, steps, εf, εx, limit, console=True):
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
                        yield bisect(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, sense=sense)
                    if method == Solver.NT:
                        yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense)
        x_prev = x.val
        f0_prev = f.val

class Context:
    places = 3

def mplot(model, x_min=-8.0, x_max=8.0, steps=1000, y_min=-10.0, y_max=10.0, nt=False):
    solver = Solver.NT if nt else Solver.BI
    for result in analyze(model, solver, x_min, x_max, steps, εf=1e-9, εx=1e-9, limit=101):
        if result.count < 101:
            print(result, file=stderr)
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Variable', color='.2')
    ax1.set_ylabel(f'Function value and derivative', color='.2')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = [[] for _ in range(3)]
    step = (x_max - x_min) / (steps - 1)
    for k in range(steps):
        x = Dual.get(x_min + k * step).var
        y = model(x)
        for d_term, p_term in zip(data, [x.val] + [y.val] + [y.der]):
            d_term.append(p_term)
    for c in reversed(range(1, 3)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')
    pyplot.show()

if __name__ == "__main__":
    mplot(lambda a: (a - 1)**2 / (a.cosh + 1).ln - 1)
else:
    print(__name__ + " module loaded", file=stderr)
