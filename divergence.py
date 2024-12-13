#!/usr/bin/env python3
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
import print_args
from math import sqrt
from sys import argv, stderr

def diverged(x, y, thresh):
    return sqrt(
        (float(x[0]) - float(y[0]))**2 + (float(x[1]) - float(y[1]))**2 + (float(x[2]) - float(y[2]))**2) > thresh

def main():
    if len(argv) < 4:
        raise Exception(">>> ERROR! Please supply two file names, a time column, and a list of thresholds <<<")
    with open(argv[1]) as a, open(argv[2]) as b:
        line_a, line_b = a.readline(), b.readline()
        for threshold in (float(thresh) for thresh in argv[3:]):
            while line_a and line_b:
                data_a, data_b = line_a.split(), line_b.split()
                if diverged(data_a, data_b, threshold):
                    print(f"threshold: {threshold:.1e}, t: {float(data_a[3]):.3f}, cpu: {float(data_a[7]):.3f}")
                    break
                line_a, line_b = a.readline(), b.readline()


main()
