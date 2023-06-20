#!/usr/bin/env python3
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv, stdin, stderr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def data_gen():
    line = stdin.readline()
    while line:
        p = line.split()
        yield float(p[3]), float(p[0]), float(p[1]), float(p[2]), p[4], p[5], p[6]
        line = stdin.readline()

def line_data():
    line_x.set_data(t_data, x_data)
    line_y.set_data(t_data, y_data)
    line_z.set_data(t_data, z_data)
    line_xm.set_data(txm_data, xm_data)
    line_ym.set_data(tym_data, ym_data)
    line_zm.set_data(tzm_data, zm_data)
    return line_x, line_y, line_z, line_xm, line_ym, line_zm,

def init():
    ax.set_ylim(minimum, maximum)
    ax.set_xlim(0, 10)
    return line_data()

def update(data):
    t, x, y, z, xm, ym, zm, = data
    t_data.append(t)
    x_data.append(x)
    y_data.append(y)
    z_data.append(z)
    if xm == 'x' or xm == 'X':
        txm_data.append(t)
        xm_data.append(x)
    if ym == 'y' or ym == 'Y':
        tym_data.append(t)
        ym_data.append(y)
    if zm == 'z' or zm == 'Z':
        tzm_data.append(t)
        zm_data.append(z)
    x_min, x_max = ax.get_xlim()
    if t >= 0.99 * x_max:
        ax.set_xlim(x_min, 1.1 * x_max)
        ax.figure.canvas.draw()
    return line_data()


print(f'Animated ODE Plotter: {argv}', file=stderr)
if len(argv) < 2:
    raise Exception('>>> ERROR! Please supply min and max <<<')
minimum, maximum = float(argv[1]), float(argv[2])
fig, ax = plt.subplots()
(line_x,), (line_xm,), = ax.plot([], [], 'g', lw=1, ms=0), ax.plot([], [], 'g', lw=0, marker='.', ms=5)
(line_y,), (line_ym,), = ax.plot([], [], 'y', lw=1, ms=0), ax.plot([], [], 'y', lw=0, marker='.', ms=5)
(line_z,), (line_zm,), = ax.plot([], [], 'c', lw=1, ms=0), ax.plot([], [], 'c', lw=0, marker='.', ms=5)
ax.grid()
t_data, x_data, y_data, z_data, txm_data, tym_data, tzm_data, xm_data, ym_data, zm_data = [], [], [], [], [], [], [], [], [], []
_ = animation.FuncAnimation(fig, update, data_gen, blit=True, interval=10, repeat=False, init_func=init)
plt.show()
