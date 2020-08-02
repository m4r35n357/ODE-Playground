#!/usr/bin/env python3

from sys import argv, stdin, stderr
from matplotlib import pyplot

def main():
    print("X-Y Plotter: {}".format(argv), file=stderr)
    if len(argv) < 2:
        raise Exception('>>> ERROR! Please supply two quantities to plot <<<')
    coordinate_a, coordinate_b = int(argv[1]), int(argv[2])
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Column ' + str(coordinate_a), color='g')
    ax1.set_ylabel('Column ' + str(coordinate_b), color='r')
    n = 0
    x, y = [], []
    while line:
        p = line.split()
        x.append(float(p[coordinate_a]))
        y.append(float(p[coordinate_b]))
        line = stdin.readline()
        n += 1
    ax1.plot(x, y, 'bo-', markersize=0)
    pyplot.show()

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
