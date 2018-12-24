#!/usr/bin/env python3

from sys import argv, stdin, stderr

from matplotlib import pyplot


def main():
    print("X-Y Plotter: {}".format(argv))
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply x/y limits <<<')
    x_max = float(argv[1])
    y_max = float(argv[2])
    line = stdin.readline()
    columns = len(line.split())
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('x', color='.2')
    ax1.set_ylabel('function and derivatives', color='.2')
    ax1.set_xlim(-x_max, x_max)
    ax1.set_ylim(-y_max, y_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = []
    for c in range(columns):
        data.append([])
    while line:
        p = line.split()
        for c in range(columns):
            try:
                data[c].append(float(p[c]))
            except ValueError:
                print(p)
        line = stdin.readline()
    for c in range(columns - 1, 0, -1):
        ax1.plot(data[0], data[c], '{}'.format(colour[c - 1]), linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')
    pyplot.show()


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
