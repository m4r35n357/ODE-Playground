#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Example: ./plotChaos-cns.py </tmp/results
#
#  Example basename='images/auto'; t=1000; while [ $t -le 10000 ]; do echo "$basename-$t"; time -p ./chaos-scan-cns thomas .01 .21 .1 $t 2>/dev/null | ./plotChaos-cns.py "$basename-$t"; t=$((t + 1000)); done

from sys import argv, stdin, stderr
from matplotlib import pyplot

def main():
    print(f'Chaos Scan Plotter: {argv}', file=stderr)
    save = False
    if len(argv) == 2:
        save = True
    coordinate_a, coordinate_b, coordinate_c = 0, 4, 7
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('Parameter', color='g')
    ax1.set_ylabel('Separation', color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Result', color='r')
    n = 0
    x, y, z = [], [], []
    while line:
        p = line.split()
        x.append(float(p[coordinate_a]))
        y.append(float(p[coordinate_b]))
        z.append(float(p[coordinate_c]))
        line = stdin.readline()
        n += 1
    ax1.plot(x, y, 'bo-', markersize=0)
    ax2.plot(x, z, 'ro-', markersize=0)
    pyplot.savefig(argv[1]) if save else pyplot.show()

main()
