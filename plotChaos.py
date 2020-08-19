#!/usr/bin/env python3

from sys import argv, stdin, stderr
from matplotlib import pyplot

def main():
    print("W-X-Y-Z Plotter: {}".format(argv), file=stderr)
    if len(argv) < 4:
        raise Exception('>>> ERROR! Please supply four quantities to plot <<<')
    coordinate_a, coordinate_b, coordinate_c, coordinate_d = int(argv[1]), int(argv[2]), int(argv[3]), int(argv[4])
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Column ' + str(coordinate_a), color='g')
    ax1.set_ylabel('Columns ' + str(coordinate_b) + ' & ' + str(coordinate_c), color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Column ' + str(coordinate_d), color='r')
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
