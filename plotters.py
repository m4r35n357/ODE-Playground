#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from matplotlib import pyplot
from ad import Series
from playground import analyze, Solver, Mode

def mplot(model, order, x_min, x_max, steps, y_min, y_max, newton=False, mode=Mode.ALL):
    #  Example: mplot(lambda x: x**3 + 3 * x**2 - 3, 5, -8, 8, 100, -10, 10)
    solver = Solver.NT if newton else Solver.BI
    for result in analyze(model, solver, x_min, x_max, steps, 1e-9, 1e-9, limit=101, order=order, mode=mode):
        if result.count < 101:
            print(result)
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('x', color='.2')
    ax1.set_ylabel('function and derivatives', color='.2')
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = []
    for c in range(order + 1):
        data.append([])
    step = (x_max - x_min) / (steps - 1)
    for k in range(steps):
        x = Series.get(order, x_min + k * step).var  # make x a variable to see derivatives!
        p = [x.val] + (~ model(x)).jet
        for c in range(order + 1):
            data[c].append(p[c])
    for c in reversed(range(order + 1)[1:]):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')
    pyplot.show()

print(__name__ + " module loaded", file=stderr)
