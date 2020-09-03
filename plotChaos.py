#!/usr/bin/env python3

from sys import argv, stdin, stderr
from matplotlib import pyplot

def main():
    print("Chaos Plotter: {}".format(argv), file=stderr)
    if len(argv) > 1:
        raise Exception('>>> ERROR! No parameters, thanks! <<<')
    coordinate_a, coordinate_b, coordinate_c, coordinate_d = 0, 4, 5, 8
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Parameter', color='g')
    ax1.set_ylabel('Separations', color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Ratio', color='r')
    n = 0
    w, x, y, z = [], [], [], []
    while line:
        p = line.split()
        w.append(float(p[coordinate_a]))
        x.append(float(p[coordinate_b]))
        y.append(float(p[coordinate_c]))
        z.append(float(p[coordinate_d]))
        line = stdin.readline()
        n += 1
    ax1.plot(w, y, 'co-', markersize=0)
    ax1.plot(w, x, 'bo-', markersize=0)
    ax2.plot(w, z, 'ro-', markersize=0)
    pyplot.show()

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
