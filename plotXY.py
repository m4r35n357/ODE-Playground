#!/usr/bin/env python3

from sys import argv, stdin, stderr

from matplotlib import pyplot


def main():
    print("X-Y Plotter: {}".format(argv))
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply a plotting interval and two quantities to plot <<<')
    interval = int(argv[1])
    coordinate1 = int(argv[2])
    coordinate2 = int(argv[3])
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Column ' + str(coordinate1), color='g')
    ax1.set_ylabel('Column ' + str(coordinate2), color='r')
    n = 0
    x = []
    y = []
    while line:
        p = line.split(' ')
        if n % interval == 0:
            x.append(float(p[coordinate1]))
            y.append(float(p[coordinate2]))
        line = stdin.readline()
        n += 1
    ax1.plot(x, y, 'bo-', markersize=0)
    pyplot.show()


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
