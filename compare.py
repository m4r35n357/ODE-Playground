#!/usr/bin/env python3

# Example: ./compare.py /tmp/dataA /tmp/dataB 3

from sys import argv, stderr
from math import sqrt
from matplotlib import pyplot

def position_error(a, b):
    return sqrt((float(a[0]) - float(b[0]))**2 + (float(a[1]) - float(b[1]))**2 + (float(a[2]) - float(b[2]))**2)

def main():
    print(f'Compare: {argv}', file=stderr)
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please supply two file names and a time coordinate column <<<")
    time = int(argv[3])
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('time', color='b')
    ax1.set_ylabel('x, y, z', color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('error', color='b')
    with open(argv[1]) as a, open(argv[2]) as b:
        line_a, line_b = a.readline(), b.readline()
        s, t, u, v, w, x, y, z = [], [], [], [], [], [], [], []
        while line_a and line_b:
            if 'nan' in line_a or 'nan' in line_b:
                break
            data_a, data_b = line_a.split(), line_b.split()
            s.append(float(data_a[time]))
            t.append(float(data_a[0]))
            u.append(float(data_b[0]))
            v.append(float(data_a[1]))
            w.append(float(data_b[1]))
            x.append(float(data_a[2]))
            y.append(float(data_b[2]))
            z.append(position_error(data_a, data_b))
            line_a, line_b = a.readline(), b.readline()
    ax1.plot(s, t, 'ko-', linewidth=1, markersize=0)
    ax1.plot(s, u, 'go-', linewidth=1, markersize=0)
    ax1.plot(s, v, 'ko-', linewidth=1, markersize=0)
    ax1.plot(s, w, 'yo-', linewidth=1, markersize=0)
    ax1.plot(s, x, 'ko-', linewidth=1, markersize=0)
    ax1.plot(s, y, 'co-', linewidth=1, markersize=0)
    ax2.plot(s, z, 'ro-', linewidth=1, markersize=0)
    pyplot.show()

main()
