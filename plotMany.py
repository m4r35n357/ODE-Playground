#!/usr/bin/env python3

from sys import argv, stdin, stderr

from matplotlib import pyplot


def main():
    print("X-Y Plotter: {}".format(argv))
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply order and x/y limits <<<')
    columns = int(argv[1])
    x_max = float(argv[2])
    y_max = float(argv[3])
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('x', color='.2')
    ax1.set_ylabel('function and derivatives', color='.2')
    ax1.set_xlim(-x_max, x_max)
    ax1.set_ylim(-y_max, y_max)
    colour = ['k', 'r', 'g', 'b', 'y', 'c', 'm']
    data = []
    for c in range(columns):
        data.append([])
    while line:
        p = line.split(' ')
        for c in range(columns):
            data[c].append(float(p[c]))
        line = stdin.readline()
    for c in range(1, columns):
        ax1.plot(data[0], data[c], '{}-'.format(colour[c - 1]), linewidth=1, markersize=0)
    pyplot.show()


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
