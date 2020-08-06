#!/usr/bin/env python3

# Example: ./ic .00001 ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .125 .1
#          ./chaos-distance.py /tmp/dataA /tmp/dataB /tmp/dataC /tmp/dataD /tmp/dataE /tmp/dataF /tmp/dataG 50001 0.000001

from sys import argv, stderr
from math import sqrt, log

RED = '\x1b[1;31m'
GREEN = '\x1b[1;32m'
ORANGE = '\x1b[0;33m'
YELLOW = '\x1b[1;33m'
BLUE = '\x1b[1;34m'
MAGENTA = '\x1b[1;35m'
CYAN = '\x1b[1;36m'
GREY = '\x1B[0m\x1b[2;37m'
WHITE = '\x1b[1;37m'
NORMAL = '\x1B[0m'

def line_to_data(line):
    data = []
    for number in line.split():
        data.append(float(number))
    return data

def worst_deviation(ref, trajectory):
    d_max = 0.0
    for point in trajectory:
        dist = sqrt((ref[0] - point[0])**2 + (ref[1] - point[1])**2 + (ref[2] - point[2])**2)
        d_max = dist if dist > d_max else d_max
    return d_max

def clip (a, limit, below=False):
    if below:
        return a if a > limit else limit
    return a if a < limit else limit

def read_files(file_a, file_b, file_c, file_d, file_e, file_f, file_g):
    data_a, data_b, data_c, data_d, data_e, data_f, data_g = [], [], [], [], [], [], []
    with open(file_a) as a, open(file_b) as b, open(file_c) as c, open(file_d) as d, open(file_e) as e, open(file_f) as f, open(file_g) as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a.append(line_to_data(line_a))
            data_b.append(line_to_data(line_b))
            data_c.append(line_to_data(line_c))
            data_d.append(line_to_data(line_d))
            data_e.append(line_to_data(line_e))
            data_f.append(line_to_data(line_f))
            data_g.append(line_to_data(line_g))
    return data_a, data_b, data_c, data_d, data_e, data_f, data_g

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the intended data file length, and two separations <<<")
    # read first set of data
    data_a1, data_b1, data_c1, data_d1, data_e1, data_f1, data_g1 = read_files('/tmp/dataA1', '/tmp/dataB1', '/tmp/dataC1', '/tmp/dataD1', '/tmp/dataE1', '/tmp/dataF1', '/tmp/dataG1')
    # read second set of data
    data_a2, data_b2, data_c2, data_d2, data_e2, data_f2, data_g2 = read_files('/tmp/dataA2', '/tmp/dataB2', '/tmp/dataC2', '/tmp/dataD2', '/tmp/dataE2', '/tmp/dataF2', '/tmp/dataG2')
    # process data sets
    length_g1, length_g2 = len(data_g1), len(data_g2)
    data_length = int(argv[1])
    separation1, separation2 = float(argv[2]), float(argv[3])
    slope = separation1 / separation2
    wd1, wd2, d1, d2 = 0.0, 0.0, 0.0, 0.0
    if length_g1 == data_length and length_g2 == data_length:
        for r, t in zip(data_g1, zip(data_a1, data_b1, data_c1, data_d1, data_e1, data_f1)):
            d1 = worst_deviation(r, t)
            wd1 = d1 if d1 > wd1 else wd1
        for r, t in zip(data_g2, zip(data_a2, data_b2, data_c2, data_d2, data_e2, data_f2)):
            d2 = worst_deviation(r, t)
            wd2 = d2 if d2 > wd2 else wd2
    # analyze data
    if length_g1 < data_length or length_g2 < data_length:
        print(f'  {WHITE}UNBOUNDED{NORMAL} value = {10.0:.1f} {10.0:.1f} ratio = {-100.0:.1f} {length_g1} / {length_g2} lines out of {data_length}')
    elif d1 < separation1 and d2 < separation2:
        d1, d2 = clip(d1, 1e-21, below=True), clip(d2, 1e-21, below=True)
        print(f'  {BLUE}CONVERGED{NORMAL} value = {log(d1):.3e} {log(d2):.3e} ratio = {clip(d1 / d2, 1500.0):.1f}')
    elif 0.8 * slope < d1 / d2 < 1.2 * slope:
        print(f'{GREEN}LIMIT-CYCLE{NORMAL} value = {log(d1):.3e} {log(d2):.3e} ratio = {clip(d1 / d2, 1500.0):.1f}')
    elif 0.5 < wd1 / wd2 < 2.0:
        print(f'    {RED}CHAOTIC{NORMAL} value = {log(wd1):.3e} {log(wd2):.3e} ratio = {clip(wd1 / wd2, 1500.0):.1f}')
    else:
        label = f"       {MAGENTA}HIGH" if wd1 / wd2 > 1000.0 else (f"        {YELLOW}MID" if 1.0 < wd1 / wd2 < 1000.0 else f"        {CYAN}LOW")
        print(f'{label}{NORMAL} value = {log(wd1):.3e} {log(wd2):.3e} ratio = {clip(wd1 / wd2, 1500.0):.1f}')

print(f'SCAN: {argv}', file=stderr)
scan()
