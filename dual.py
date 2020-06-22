#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from math import sqrt, sin, cos, sinh, cosh, tan, tanh, exp, log, copysign
from matplotlib import pyplot

class Context:
    places = 3

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
        return self + o

    def __sub__(self, o):
        return self + (- o)

    def __rsub__(self, o):
        return - self + o

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val * o.val, self.der * o.val + self.val * o.der)
        elif isinstance(o, (float, int)):
            return Dual(self.val * o, self.der * o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self * o

    def __truediv__(self, o):
        if isinstance(o, Dual):
            assert abs(o.val) != 0.0, f"other.val = {o.val}"
            return Dual(self.val / o.val, (self.der * o.val - self.val * o.der) / o.val**2)
        elif isinstance(o, (float, int)):
            assert abs(o) != 0.0, f"other = {o}"
            return Dual(self.val / o, self.der / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        assert self.val != 0.0, f"self.val = {self.val}"
        return Dual(o / self.val, - self.der * o / self.val**2)

    def __pow__(self, o):
        if isinstance(o, int):
            i_pow = Dual(self.val, self.der)
            for _ in range(abs(o) - 1):
                i_pow = i_pow * self
            return i_pow if o > 0 else (1.0 / i_pow if o < 0 else Dual(1.0, 0.0))
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
            return Dual(self.val**o, self.der * o * self.val**(o - 1))
        elif isinstance(o, Dual):
            return (self.ln * o).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self * log(o)).exp

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
        t = tan(self.val)
        return Dual(t, self.der * (1.0 + t**2))

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.der * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.der * sinh(self.val))

    @property
    def tanh(self):
        t = tanh(self.val)
        return Dual(t, self.der * (1.0 - t**2))

    @property
    def sqr(self):
        return Dual(self.val * self.val, self.der * 2.0 * self.val)

    @property
    def sqrt(self):
        assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
        r = sqrt(self.val)
        return Dual(r, self.der * 0.5 / r)

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

class Matrix3x3(namedtuple('Matrix3x3Type', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])):
    def __str__(self):
        return f'{self.a:+}, {self.b:+}, {self.c:+}\n{self.d:+}, {self.e:+}, {self.f:+}\n{self.g:+}, {self.h:+}, {self.i:+}'

def function_3(model_a, model_b, model_c, x, y, z):
    a, b, c = Dual.get(x), Dual.get(y), Dual.get(z)
    return model_a(a, b, c).val, model_b(a, b, c).val, model_c(a, b, c).val

def gradient_3(model, x, y, z):
    a, b, c = Dual.get(x), Dual.get(y), Dual.get(z)
    return model(a.var, b, c).der, model(a, b.var, c).der, model(a, b, c.var).der

def jacobian(model_a, model_b, model_c, x, y, z):
    a, b, c = Dual.get(x), Dual.get(y), Dual.get(z)
    return Matrix3x3(a=model_a(a.var, b, c).der, b=model_a(a, b.var, c).der, c=model_a(a, b, c.var).der,
                     d=model_b(a.var, b, c).der, e=model_b(a, b.var, c).der, f=model_b(a, b, c.var).der,
                     g=model_c(a.var, b, c).der, h=model_c(a, b.var, c).der, i=model_c(a, b, c.var).der)

def determinant_3x3(m):
    c = Matrix3x3(a=(m.e * m.i - m.f * m.h), b=-(m.d * m.i - m.f * m.g), c=(m.d * m.h - m.e * m.g),
                  d=-(m.b * m.i - m.c * m.h), e=(m.a * m.i - m.c * m.g), f=-(m.a * m.h - m.b * m.g),
                  g=(m.b * m.f - m.c * m.e), h=-(m.a * m.f - m.c * m.d), i=(m.a * m.e - m.b * m.d))
    d =  m.a * c.a + m.b * c.b + m.c * c.c
    assert d != 0.0, f"determinant = {d}"
    return d, c

def invert_3x3(m, transpose=True):
    d, c = determinant_3x3(m)
    if transpose:
        return Matrix3x3(a=c.a/d, b=c.d/d, c=c.g/d,
                         d=c.b/d, e=c.e/d, f=c.h/d,
                         g=c.c/d, h=c.f/d, i=c.i/d)
    return Matrix3x3(a=c.a/d, b=c.b/d, c=c.c/d,
                     d=c.d/d, e=c.e/d, f=c.f/d,
                     g=c.g/d, h=c.h/d, i=c.i/d)

def c_e(m, l):
    lam = Dual.get(l).var
    return determinant_3x3(Matrix3x3(a=m.a - lam, b=m.b, c=m.c,
                                     d=m.d, e=m.e - lam, f=m.f,
                                     g=m.g, h=m.h, i=m.i - lam))[0]

def equilibrium(model_a, model_b, model_c, x, y, z, max_it=100):
    fx, fy, fz = function_3(model_a, model_b, model_c, x, y, z)
    count = 0
    print(count, x, y, z, fx, fy, fz)
    while abs(fx) + abs(fy) + abs(fz) > 1e-12 and count < max_it:
        j_1 = invert_3x3(jacobian(model_a, model_b, model_c, x, y, z))
        x -= j_1.a * fx + j_1.b * fy + j_1.c * fz
        y -= j_1.d * fx + j_1.e * fy + j_1.f * fz
        z -= j_1.g * fx + j_1.h * fy + j_1.i * fz
        fx, fy, fz = function_3(model_a, model_b, model_c, x, y, z)
        count += 1
        print(count, x, y, z, fx, fy, fz)
    return x, y, z, fx, fy, fz

def ode_model(l):
    return l

def analyze_3(model, model_a, model_b, model_c, x, y, z, x0, x1, steps, εf, εx, limit, console=True):
    f0_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        l = Dual.get(x0 + k * step).var
        j = jacobian(model_a, model_b, model_c, x, y, z)
        f = determinant_3x3(Matrix3x3(a=j.a - l, b=j.b, c=j.c,
                                     d=j.d, e=j.e - l, f=j.f,
                                     g=j.g, h=j.h, i=j.i - l))[0]
        if not console:
            print(f'{x.val:.{Context.places}e} {f}')
        if k > 0:
            if f0_prev * f.val < 0.0:
                sense = Sense.DECREASING if f0_prev > f.val else Sense.INCREASING
                yield newton(model, x.val, εf=εf, εx=εx, limit=limit, sense=sense)
        f0_prev = f.val

def bisect(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    a, b, c = Dual.get(xa), Dual.get(xb), Dual.get(3)
    f_sign = copysign(1, model(Dual.get(xa)).val)
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
            print(Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=fc, δx=δx), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=fc, δx=δx)

def newton(model, x0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    x, f = Dual.get(x0).var, Dual.get(1)
    δx = i = 1
    while i <= limit and abs(f.val) > εf or abs(δx) > εx:
        f = model(x)
        δx = - f.val / f.der
        x += δx
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=x.val, f=f.val, δx=δx), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=x.val, f=f.val, δx=δx)

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

def mplot(model, x_min=-8.0, x_max=8.0, steps=1000, y_min=-10.0, y_max=10.0, nt=False):
    solver = Solver.NT if nt else Solver.BI
    for result in analyze(model, solver, x_min, x_max, steps, εf=1e-9, εx=1e-9, limit=101):
        # noinspection PyTypeChecker
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
