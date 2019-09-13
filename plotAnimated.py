#!/usr/bin/env python3

from sys import argv, stdin, stderr
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def data_gen():
    line = stdin.readline()
    while line:
        p = line.split(' ')
        yield float(p[3]), float(p[0]), float(p[1]), float(p[2])
        line = stdin.readline()


def line_data():
    line_x.set_data(t_data, x_data)
    line_y.set_data(t_data, y_data)
    line_z.set_data(t_data, z_data)
    return line_x, line_y, line_z,


def init():
    ax.set_ylim(minimum, maximum)
    ax.set_xlim(0, 10)
    return line_data()


def update(data):
    t, x, y, z = data
    t_data.append(t)
    x_data.append(x)
    y_data.append(y)
    z_data.append(z)
    x_min, x_max = ax.get_xlim()
    if t >= 0.99 * x_max:
        ax.set_xlim(x_min, 1.1 * x_max)
        ax.figure.canvas.draw()
    return line_data()


interval, minimum, maximum = int(argv[1]), float(argv[2]), float(argv[3])
fig, ax = plt.subplots()
line_x, = ax.plot([], [], 'g', lw=1)
line_y, = ax.plot([], [], 'y', lw=1)
line_z, = ax.plot([], [], 'c', lw=1)
ax.grid()
t_data, x_data, y_data, z_data = [], [], [], []


if __name__ == "__main__":
    _ = animation.FuncAnimation(fig, update, data_gen, blit=True, interval=10, repeat=False, init_func=init)
    plt.show()
else:
    print(__name__ + " module loaded", file=stderr)
