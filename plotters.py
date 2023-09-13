#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from collections import namedtuple
from enum import Enum, unique
from matplotlib import pyplot
from ad import Series, Dual, Context


@unique
class Solver(Enum):
    NA = "No analysis, or all"
    BI = "Bisection method"
    NT = "Newton-Raphson method"

@unique
class Label(Enum):
    MIN = '/'
    MAX = '\\'
    FLAT = '_'

@unique
class Mode(Enum):
    ROOT___ = 0
    MIN_MAX = 1
    INFLECT = 2
    ALL = 3


class Result(namedtuple('ResultType', ['method', 'x', 'f', 'δx', 'count', 'label', 'mode'])):
    def __str__(self):
        return f'{self.method}  x: {self.x:+.{Context.places}e}  δx: {self.δx:+.{Context.places}e}  ' \
               f'f: {self.f:+.{Context.places}e}  {self.label} {self.mode} {self.count}'

def _analyze_s(model, method, x0, x1, steps, εf, εx, limit, mode, console, debug):
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
                    sense = Label.MAX if f0_prev > f.val else Label.MIN
                    if method == Solver.BI:
                        yield bisect_s(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, label=sense, debug=debug)
                    if method == Solver.NT:
                        yield newton_s(model, x.val, εf=εf, εx=εx, limit=limit, label=sense, debug=debug)
                if (mode == Mode.MIN_MAX or mode == Mode.ALL) and f1_prev * f.jet[Mode.MIN_MAX.value] < 0.0:
                    sense = Label.MAX if f1_prev > f.jet[Mode.MIN_MAX.value] else Label.MIN
                    if method == Solver.BI:
                        yield bisect_s(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, label=sense, mode=Mode.MIN_MAX, debug=debug)
                    if method == Solver.NT:
                        yield newton_s(model, x.val, εf=εf, εx=εx, limit=limit, label=sense, mode=Mode.MIN_MAX, debug=debug)
                if (mode == Mode.INFLECT or mode == Mode.ALL) and f2_prev * f.jet[Mode.INFLECT.value] < 0.0:
                    sense = Label.MAX if f2_prev > f.jet[Mode.INFLECT.value] else Label.MIN
                    if method == Solver.BI:
                        yield bisect_s(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, label=sense, mode=Mode.INFLECT, debug=debug)
                    if method == Solver.NT:
                        yield newton_s(model, x.val, εf=εf, εx=εx, limit=limit, label=sense, mode=Mode.INFLECT, debug=debug)
        x_prev = x.val
        f0_prev = f.val
        f1_prev = f.jet[Mode.MIN_MAX.value]
        f2_prev = f.jet[Mode.INFLECT.value]

def _analyze_d(model, method, x0, x1, steps, εf, εx, limit, console, debug):
    x_prev = f0_prev = None
    step = (x1 - x0) / (steps - 1)
    for k in range(steps):
        x = Dual(x0 + k * step).var
        f = model(x)
        if not console:
            print(f'{x.val:.{Context.places}e} {f}')
        if method != Solver.NA:
            if k > 0:
                if f0_prev * f.val < 0.0:
                    sense = Label.MAX if f0_prev > f.val else Label.MIN
                    if method == Solver.BI:
                        yield bisect_d(model, x.val, x_prev, εf=εf, εx=εx, limit=limit, label=sense, debug=debug)
                    if method == Solver.NT:
                        yield newton_d(model, x.val, εf=εf, εx=εx, limit=limit, label=sense, debug=debug)
        x_prev = x.val
        f0_prev = f.val

def _plot_s(model, order, x_min, x_max, steps, y_min, y_max):
    #  Plot the function and its derivatives
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('x.val', color='.2')
    ax1.set_ylabel(f'f(x) and the first {order} derivatives', color='.2')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = [[] for _ in range(order + 2)]
    step = (x_max - x_min) / (steps - 1)
    for k in range(1, steps - 1):
        x = Series.get(order + 1, x_min + k * step)
        x = x.var if order > 0 else x
        for d_term, p_term in zip(data, [x.val] + (~ model(x)).jet):
            d_term.append(p_term)
    for c in reversed(range(1, order + 2)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')

def _plot_d(model, x_min, x_max, steps, y_min, y_max):
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Variable', color='.2')
    ax1.set_ylabel(f'Function value and derivative', color='.2')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = [[] for _ in range(3)]
    step = (x_max - x_min) / (steps - 1)
    for k in range(1, steps - 1):
        x = Dual(x_min + k * step).var
        y = model(x)
        for d_term, p_term in zip(data, [x.val] + [y.val] + [y.dot]):
            d_term.append(p_term)
    for c in reversed(range(1, 3)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')

def bisect_s(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, label=Label.FLAT, mode=Mode.ROOT___, debug=False):
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
            b, a = c, c
        δx = b.val - a.val
        i += 1
        if debug:
            print(Result(method=Solver.BI.name, count=i-1, label=label.value, mode=mode.name, x=c.val, f=fc, δx=δx), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, label=label.value, mode=mode.name, x=c.val, f=fc, δx=δx)

def bisect_d(model, xa, xb, εf=1e-12, εx=1e-12, limit=101, label=Label.FLAT, debug=False):
    a, b, c = Dual(xa), Dual(xb), Dual(3)
    f_sign = model(Dual(xa)).val
    δx = fc = i = 1
    while i <= limit and (abs(fc) > εf or abs(δx) > εx):
        c = 0.5 * (a + b)
        fc = model(c).val
        if f_sign * fc < 0.0:
            b = c
        elif f_sign * fc > 0.0:
            a = c
        else:
            b, a = c, c
        δx = b.val - a.val
        i += 1
        if debug:
            print(Result(method=Solver.BI.name, count=i-1, label=label.value, x=c.val, f=fc, δx=δx, mode=Mode.ROOT___.name), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, label=label.value, x=c.val, f=fc, δx=δx, mode=Mode.ROOT___.name)

def newton_s(model, x0, εf=1e-12, εx=1e-12, limit=101, label=Label.FLAT, mode=Mode.ROOT___, debug=False):
    x, f = Series.get(2 + mode.value, x0).var, [1.0, 0.0]
    δx = i = 1
    while i <= limit and (abs(f[0]) > εf or abs(δx) > εx):
        f = (~ model(x)).jet[mode.value:2 + mode.value]
        δx = - f[0] / f[1]
        x += δx
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, label=label.value, mode=mode.name, x=x.val, f=f[0], δx=δx), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, label=label.value, mode=mode.name, x=x.val, f=f[0], δx=δx)

def newton_d(model, x0, εf=1e-12, εx=1e-12, limit=101, label=Label.FLAT, debug=False):
    x, f = Dual(x0).var, Dual(1)
    δx = i = 1
    while i <= limit and (abs(f.val) > εf or abs(δx) > εx):
        f = model(x)
        δx = - f.val / f.dot
        x += δx
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, label=label.value, x=x.val, f=f.val, δx=δx, mode=Mode.ROOT___.name), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, label=label.value, x=x.val, f=f.val, δx=δx, mode=Mode.ROOT___.name)

def scan_s(model, x_min=-2.0, x_max=2.0, steps=1000, εf=1e-9, εx=1e-9, limit=101, newton=True, mode=Mode.ALL, console=True, debug=False):
    for result in _analyze_s(model, Solver.NT if newton else Solver.BI, x_min, x_max, steps, εf, εx, limit, mode, console, debug):
        if result.count < 101:
            print(result, file=stderr)

def scan_d(model, x_min=-2.0, x_max=2.0, steps=1000, εf=1e-9, εx=1e-9, limit=101, newton=True, console=True, debug=False):
    for result in _analyze_d(model, Solver.NT if newton else Solver.BI, x_min, x_max, steps, εf, εx, limit, console, debug):
        if result.count < 101:
            print(result, file=stderr)

def mplot_s(model, order=12, x_min=-2.0, x_max=2.0, steps=1000, y_min=-10.0, y_max=10.0):
    _plot_s(model, order, x_min, x_max, steps, y_min, y_max)
    pyplot.show()

def mplot_d(model, x_min=-2.0, x_max=2.0, steps=1000, y_min=-10.0, y_max=10.0):
    _plot_d(model, x_min, x_max, steps, y_min, y_max)
    pyplot.show()

def msave_s(filename, model, order=12, x_min=-2.0, x_max=2.0, steps=1000, y_min=-10.0, y_max=10.0):
    _plot_s(model, order, x_min, x_max, steps, y_min, y_max)
    pyplot.savefig(filename)

def msave_d(filename, model, x_min=-2.0, x_max=2.0, steps=1000, y_min=-10.0, y_max=10.0):
    _plot_d(model, x_min, x_max, steps, y_min, y_max)
    pyplot.savefig(filename)


if __name__ == "__main__":  # hard-coded example
    function = lambda a: (a - 1)**2 / (a.cosh + 1).ln - 1
    print(f'Root Scan (Dual)')
    scan_d(function)
    mplot_d(function)
    print(f'Multi Scan (Series)')
    scan_s(function)
    mplot_s(function)
else:
    print(__name__ + " module loaded", file=stderr)

'''
from ad import *
from plotters import *

a = Series.get(4, 3.0).var
print(a)
print(~(a * a))

b = Dual(3.0).var
print(b)
print(b * b)

f = lambda a: a * a * a + 2.0 * a * a - 3.0 * a + 1.0

scan_s(f)
mplot_s(f)

scan_d(f)
mplot_d(f)


bisect_d(f, -4.0, -2.0)
newton_d(f, -4.0)

bisect_s(f, -4.0, -2.0)
newton_s(f, -4.0)


f = lambda a: (a.exp + (a.sqr - 4.0).exp).ln - value
f = lambda a: (a.sqr + (a.exp - 4).sqr).sqrt - value
f = lambda a: (2 * a).sin - 2 * a.sin * a.cos - value
f = lambda a: (2 * a).cos - a.cos.sqr + a.sin.sqr - value
f = lambda a: (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
f = lambda a: (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
f = lambda a: a.sqr.sqrt - value
f = lambda a: a.exp.ln - value
f = lambda a: a.tan - a.sin / a.cos - value
f = lambda a: a.tanh - a.sinh / a.cosh - value
f = lambda a: a.sin.asin - value
f = lambda a: a.cos.acos - value
f = lambda a: a.tan.atan - value
f = lambda a: (a + 7) / (3.0 - a)

for i in range(-50, 51):
    print(i * 0.1)
    f = lambda a: (a.sqr + 1)**(i*0.1)
    #scan(f)
    #mplot(f)
    msave_s(f'output/test_{i+50:03d}.png', f)

for i in range(19):
    f = lambda a: (a.sqr + 1 / (10**i)).sqrt
    msave_s(f'output/test_{i:03d}.png', f)

f = lambda a: abs(- (a - 1.5)**3 + 4 * (a - 1.5) - 2)
f = lambda a: - (abs(a) - 1.5)**3 + 4 * (abs(a) - 1.5) - 2

f = lambda a: a.sin / pi + (3 * a).sin / (3 * pi) + (5 * a).sin / (5 * pi) + (7 * a).sin / (7 * pi)

rm -f output/*.png
ffmpeg -y -i output/test_%03d.png taylor.mp4
mplayer -fps 1 taylor.mp4
'''
