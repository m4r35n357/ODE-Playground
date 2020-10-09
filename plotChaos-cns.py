#!/usr/bin/env python3

from sys import argv, stdin, stderr
from matplotlib import pyplot

def main():
    print("Chaos Plotter: {}".format(argv), file=stderr)
    if len(argv) > 1:
        raise Exception('>>> ERROR! No parameters, thanks! <<<')
    coordinate_a, coordinate_b, coordinate_c = 0, 1, 4
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Parameter', color='g')
    ax1.set_ylabel('Separation', color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Result', color='r')
    n = 0
    w, x, y, z = [], [], [], []
    while line:
        p = line.split()
        x.append(float(p[coordinate_a]))
        y.append(float(p[coordinate_c]))
        z.append(1.0 if 'CHAOTIC' in str(p[coordinate_b]) else (0.5 if 'UNCLASSIFIED' in str(p[coordinate_b]) else (-0.1 if 'UNBOUNDED' in str(p[coordinate_b]) else 0.0)))
        line = stdin.readline()
        n += 1
    ax1.plot(x, y, 'bo-', markersize=0)
    ax2.plot(x, z, 'ro-', markersize=0)
    pyplot.show()

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
