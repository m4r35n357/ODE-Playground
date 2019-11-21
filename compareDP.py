#!/usr/bin/env python3

# Example: ./compare.py /tmp/dataA /tmp/dataB 3

from math import sqrt
from sys import argv, stderr
from matplotlib import pyplot

def position_error(a, b):
    return sqrt((float(a[0]) - float(b[0]))**2 + (float(a[1]) - float(b[1]))**2)\
           + sqrt((float(a[2]) - float(b[2]))**2 + (float(a[3]) - float(b[3]))**2)

def main():
    print(f'Compare: {argv}')
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please supply two file names and a time coordinate column <<<")
    time = int(argv[3])
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('time', color='b')
    ax1.set_ylabel('x, y', color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('z, error', color='b')
    with open(argv[1]) as a, open(argv[2]) as b:
        line_a = a.readline()
        line_b = b.readline()
        t, a0, b0, a1, b1, a2, b2, a3, b3, error = [], [], [], [], [], [], [], [], [], []
        while line_a and line_b:
            if 'nan' in line_a or 'nan' in line_b:
                break
            data_a = line_a.split(' ')
            data_b = line_b.split(' ')
            t.append(float(data_a[time]))
            a0.append(float(data_a[0]))
            b0.append(float(data_b[0]))
            a1.append(float(data_a[1]))
            b1.append(float(data_b[1]))
            a2.append(float(data_a[2]))
            b2.append(float(data_b[2]))
            a3.append(float(data_a[3]))
            b3.append(float(data_b[3]))
            error.append(position_error(data_a, data_b))
            line_a = a.readline()
            line_b = b.readline()
    ax1.plot(t, a0, 'ko-', linewidth=1, markersize=0)
    ax1.plot(t, b0, 'go-', linewidth=1, markersize=0)
    ax1.plot(t, a1, 'ko-', linewidth=1, markersize=0)
    ax1.plot(t, b1, 'yo-', linewidth=1, markersize=0)
    ax2.plot(t, a2, 'ko-', linewidth=1, markersize=0)
    ax2.plot(t, b2, 'co-', linewidth=1, markersize=0)
    ax2.plot(t, a3, 'ko-', linewidth=1, markersize=0)
    ax2.plot(t, b3, 'mo-', linewidth=1, markersize=0)
    ax2.plot(t, error, 'ro-', linewidth=1, markersize=0)
    pyplot.show()

main()
