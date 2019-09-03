#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from matplotlib import pyplot
from ad import Series
from playground import analyze, Solver, Mode

def mplot(model, order=13, x_min=-8.0, x_max=8.0, steps=1000, y_min=-10.0, y_max=10.0, newton=False, mode=Mode.ALL):
    #  Find roots, turning points and inflections of the model function
    solver = Solver.NT if newton else Solver.BI
    for result in analyze(model, solver, x_min, x_max, steps, 1e-9, 1e-9, limit=101, order=order, mode=mode):
        if result.count < 101:
            print(result, file=stderr)
    #  Plot the function and its derivatives
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Variable', color='.2')
    ax1.set_ylabel(f'Function value and first {order - 1} derivatives', color='.2')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = [[] for _ in range(order + 1)]
    step = (x_max - x_min) / (steps - 1)
    for k in range(steps):
        x = Series.get(order, x_min + k * step).var
        for d_term, p_term in zip(data, [x.val] + (~ model(x)).jet):
            d_term.append(p_term)
    for c in reversed(range(1, order + 1)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')
    pyplot.show()

if __name__ == "__main__":
    mplot(lambda a: (a - 1)**2 / (a.cosh + 1).ln - 1, order=4)
else:
    print(__name__ + " module loaded", file=stderr)
