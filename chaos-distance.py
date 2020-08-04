#!/usr/bin/env python3

# Example: ./ic .00001 ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .125 .1
#          ./chaos-distance.py /tmp/dataA /tmp/dataB /tmp/dataC /tmp/dataD /tmp/dataE /tmp/dataF /tmp/dataG 50001 0.000001

from sys import argv, stderr
from math import sqrt

RED = '\x1b[1;31m'
GREEN = '\x1b[0;32m'
ORANGE = '\x1b[0;33m'
YELLOW = '\x1b[1;33m'
BLUE = '\x1b[1;34m'
CYAN = '\x1b[0;36m'
GREY = '\x1B[0m\x1b[2;37m'
WHITE = '\x1b[0;37m'
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
    data_g1_length = len(data_g1)
    # read second set of data
    data_a2, data_b2, data_c2, data_d2, data_e2, data_f2, data_g2 = read_files('/tmp/dataA2', '/tmp/dataB2', '/tmp/dataC2', '/tmp/dataD2', '/tmp/dataE2', '/tmp/dataF2', '/tmp/dataG2')
    data_g2_length = len(data_g2)
    # process data sets
    data_length = int(argv[1])
    separation1, separation2 = float(argv[2]), float(argv[3])
    slope = separation1 / separation2
    wd1, wd2 = 0.0, 0.0
    if data_g1_length == data_length and data_g2_length == data_length:
        for r, t in zip(data_g1, zip(data_a1, data_b1, data_c1, data_d1, data_e1, data_f1)):
            wd1 = worst_deviation(r, t)
        for r, t in zip(data_g2, zip(data_a2, data_b2, data_c2, data_d2, data_e2, data_f2)):
            wd2 = worst_deviation(r, t)
    # analyze data
    if data_g1_length < data_length or data_g2_length < data_length:
        print(f'  {WHITE}UNBOUNDED{NORMAL} {"dummy values"} = {-0.1:.1f} {data_g1_length} / {data_g2_length} lines out of {data_length}')
    elif wd1 < separation1 and wd2 < separation2:
        print(f'  {BLUE}CONVERGED{NORMAL} {"final values"} = {wd1:.3e} < {separation1:.1e}, {wd2:.3e} < {separation2:.1e}')
    elif 0.8 * slope < wd1 / wd2 < 1.2 * slope:
        print(f'{GREEN}LIMIT CYCLE{NORMAL} {"final values"} = {wd1:.3e} {wd2:.3e} ratio = {wd1 / wd2:.1f}')
    elif 0.5 < wd1 / wd2 < 2.0:
        print(f'    {RED}CHAOTIC{NORMAL} {"final values"} = {wd1:.3e} {wd2:.3e} ratio = {wd1 / wd2:.1f}')
    else:
        label = "       HIGH" if wd1 / wd2 > 1000.0 else ("        LOW" if 1.0 < wd1 / wd2 < 1000.0 else "   VERY LOW")
        print(f'{YELLOW}{label}{NORMAL} {"final values"} = {wd1:.3e} {wd2:.3e} ratio = {wd1 / wd2:.1f}')

print(f'SCAN: {argv}', file=stderr)
scan()
