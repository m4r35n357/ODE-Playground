#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from matplotlib import pyplot
from ad import Series, Dual, analyze_s, analyze_d, Solver, Mode


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
    for k in range(steps):
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
    for k in range(steps):
        x = Dual.get(x_min + k * step).var
        y = model(x)
        for d_term, p_term in zip(data, [x.val] + [y.val] + [y.dot]):
            d_term.append(p_term)
    for c in reversed(range(1, 3)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')

def scan(model, variable, x_min=-8.0, x_max=8.0, steps=1000, εf=1e-9, εx=1e-9, limit=101, newton=True, mode=Mode.ALL, console=True, debug=False):
    if isinstance(variable, Series):
        for result in analyze_s(model, Solver.NT if newton else Solver.BI, x_min, x_max, steps, εf, εx, limit, mode, console, debug):
            if result.count < 101:
                print(result, file=stderr)
    else:
        for result in analyze_d(model, Solver.NT if newton else Solver.BI, x_min, x_max, steps, εf, εx, limit, console, debug):
            if result.count < 101:
                print(result, file=stderr)


def mplot(model, variable, order=12, x_min=-8.0, x_max=8.0, steps=1000, y_min=-10.0, y_max=10.0):
    if isinstance(variable, Series):
        _plot_s(model, order, x_min, x_max, steps, y_min, y_max)
    else:
        _plot_d(model, x_min, x_max, steps, y_min, y_max)
    pyplot.show()

def msave(filename, model, variable, order=12, x_min=-8.0, x_max=8.0, steps=1000, y_min=-10.0, y_max=10.0):
    if isinstance(variable, Series):
        _plot_s(model, order, x_min, x_max, steps, y_min, y_max)
    else:
        _plot_d(model, x_min, x_max, steps, y_min, y_max)
    pyplot.savefig(filename)

if __name__ == "__main__":  # hard-coded example
    f = lambda a: (a - 1)**2 / (a.cosh + 1).ln - 1
    print(f'Multi Scan (Series)')
    scan_s(f)
    mplot_s(f)
    print(f'Root Scan (Dual)')
    mplot_d(f)
else:
    print(__name__ + " module loaded", file=stderr)

'''
from ad import *
from plotters import *

f = lambda a: a * a * a + 2.0 * a * a - 3.0 * a +1.0

x = Series.get(5, 1.0).var
scan(f, x)
mplot(f, x)

y = Dual.get(1.0).var
scan(f, y)
mplot(f, y)

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
    f = lambda a: (a**2 + 1)**(i*0.1)
    #scan(f)
    #mplot(f)
    msave(f'test_{i+50:03d}.png', f)

for i in range(19):
    f = lambda a: (a**2 + 1.0 / (10.0)**i)**0.5
    msave(f'test_{i:03d}.png', f)

f = lambda a: abs(- (a - 1.5)**3 + 4 * (a - 1.5) - 2)
f = lambda a: - (abs(a) - 1.5)**3 + 4 * (abs(a) - 1.5) - 2

f = lambda a: a.sin / pi + (3 * a).sin / (3 * pi) + (5 * a).sin / (5 * pi) + (7 * a).sin / (7 * pi)

ffmpeg -y -i test_%03d.png taylor.mp4
mplayer -fps 2 taylor.mp4
'''
