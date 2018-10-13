#!/usr/bin/env python3
"""
=====
Decay
=====

This example showcases a sinusoidal decay animation.
"""

from sys import argv, stdin, stderr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def data_gen(t=0):
    cnt = 0
    while cnt < 1000:
        cnt += 1
        t += 0.1
        n = 0
        data = stdin.readline()
        while data:
            p = data.split(' ')
            # if n % interval == 0:
            yield float(p[3]), float(p[0]), float(p[1]), float(p[2])
            data = stdin.readline()
            n += 1


def init():
    ax.set_ylim(minimum, maximum)
    ax.set_xlim(0, 10)
    # ax.set_autoscale_on(True)
    line_x.set_data(tdata, xdata)
    line_y.set_data(tdata, ydata)
    line_z.set_data(tdata, zdata)
    return line_x, line_y, line_z,


interval = int(argv[1])
minimum = int(argv[2])
maximum = int(argv[3])
fig, ax = plt.subplots()
line_x, = ax.plot([], [], 'g', lw=1)
line_y, = ax.plot([], [], 'y', lw=1)
line_z, = ax.plot([], [], 'c', lw=1)
ax.grid()
tdata, xdata, ydata, zdata = [], [], [], []


def run(data):
    # update the data
    t, x, y, z = data
    tdata.append(t)
    xdata.append(x)
    ydata.append(y)
    zdata.append(z)
    xmin, xmax = ax.get_xlim()

    if t >= 0.99 * xmax:
        ax.set_xlim(xmin, 1.1 * xmax)
        ax.figure.canvas.draw()
    line_x.set_data(tdata, xdata)
    line_y.set_data(tdata, ydata)
    line_z.set_data(tdata, zdata)

    return line_x, line_y, line_z,


def main():
    ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=10, repeat=False, init_func=init)
    # ax.relim()
    # ax.autoscale_view(True, True, True)
    plt.show()


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
