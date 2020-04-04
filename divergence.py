#!/usr/bin/env python3

from math import sqrt
from sys import argv, stderr

RED = '\x1b[1;31m'
GREEN = '\x1b[0;32m'
ORANGE = '\x1b[0;33m'
YELLOW = '\x1b[1;33m'
BLUE = '\x1b[1;34m'
CYAN = '\x1b[0;36m'
GREY = '\x1B[0m\x1b[2;37m'
WHITE = '\x1b[0;37m'
NORMAL = '\x1B[0m'

def diverged(x, y, thresh):
    return sqrt(
        (float(x[0]) - float(y[0]))**2 + (float(x[1]) - float(y[1]))**2 + (float(x[2]) - float(y[2]))**2) > thresh

def main():
    print(f'Divergence: {argv}')
    if len(argv) < 5:
        raise Exception(">>> ERROR! Please supply two file names, a time column, and a list of thresholds <<<")
    time = int(argv[3])
    with open(argv[1]) as a, open(argv[2]) as b:
        line_a, line_b = a.readline(), b.readline()
        for threshold in (float(thresh) for thresh in argv[4:]):
            if threshold <= 1.0e-18:
                colour = BLUE
            elif threshold <= 1.0e-12:
                colour = GREEN
            elif threshold <= 1.0e-9:
                colour = CYAN
            elif threshold <= 1.0e-6:
                colour = YELLOW
            elif threshold <= 1.0e-3:
                colour = ORANGE
            else:
                colour = RED
            while line_a and line_b:
                data_a, data_b = line_a.split(' '), line_b.split(' ')
                if diverged(data_a, data_b, threshold):
                    print(f"{GREY}Threshold: {colour}{threshold:.1e}{GREY}, t: {WHITE}{float(data_a[time]):.3f}{NORMAL}", file=stderr)
                    break
                line_a, line_b = a.readline(), b.readline()

main()
