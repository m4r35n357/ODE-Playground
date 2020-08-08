#!/usr/bin/env python3

# Example: ./ic .00001 ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .125 .1
#          ./chaos-distance.py /tmp/dataA /tmp/dataB /tmp/dataC /tmp/dataD /tmp/dataE /tmp/dataF /tmp/dataG 50001 0.000001

from sys import argv, stderr
from math import sqrt

RED = '\x1b[1;31m'
GREEN = '\x1b[1;32m'
YELLOW = '\x1b[1;33m'
BLUE = '\x1b[1;34m'
WHITE = '\x1b[1;37m'
NORMAL = '\x1B[0m'

def data(line):
    d = []
    for number in line.split():
        d.append(float(number))
    return d

def worst_separation(ref, trajectories):
    separation = 0.0
    for trajectory in trajectories:
        d = sqrt((ref[0] - trajectory[0])**2 + (ref[1] - trajectory[1])**2 + (ref[2] - trajectory[2])**2)
        separation = d if d > separation else separation
    return separation

def clip (a, limit):
    return a if a < limit else limit

def read_files(file_a, file_b, file_c, file_d, file_e, file_f, file_g):
    data_length = 0
    with open(file_a) as a, open(file_b) as b, open(file_c) as c, open(file_d) as d, open(file_e) as e, open(file_f) as f, open(file_g) as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_length += 1
    return data_length, data(line_g), (data(line_a), data(line_b), data(line_c), data(line_d), data(line_e), data(line_f))

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the expected data file length, and two initial separations <<<")
    # read in data
    length_g1, data_ref_1, data_1 = read_files('/tmp/dataA1', '/tmp/dataB1', '/tmp/dataC1', '/tmp/dataD1', '/tmp/dataE1', '/tmp/dataF1', '/tmp/dataG1')
    length_g2, data_ref_2, data_2 = read_files('/tmp/dataA2', '/tmp/dataB2', '/tmp/dataC2', '/tmp/dataD2', '/tmp/dataE2', '/tmp/dataF2', '/tmp/dataG2')
    # analyze data
    data_length, separation1, separation2 = int(argv[1]), float(argv[2]), float(argv[3])
    slope = separation1 / separation2
    if length_g1 == data_length and length_g2 == data_length:
        d1 = worst_separation(data_ref_1, data_1)
        d2 = worst_separation(data_ref_2, data_2)
        if d1 < separation1 and d2 < separation2:
            print(f'   {BLUE}CONVERGED{NORMAL} value = {d1:.3e} {d2:.3e} ratio = {slope:.1f}')
        elif 0.8 * slope < d1 / d2 < 1.2 * slope:
            print(f' {GREEN}LIMIT-CYCLE{NORMAL} value = {d1:.3e} {d2:.3e} ratio = {clip(d1 / d2, 1500.0):.1f}')
        elif 0.5 < d1 / d2 < 2.0:
            print(f'     {RED}CHAOTIC{NORMAL} value = {d1:.3e} {d2:.3e} ratio = {clip(d1 / d2, 1500.0):.1f}')
        else:
            print(f'{YELLOW}UNCLASSIFIED{NORMAL} value = {d1:.3e} {d2:.3e} ratio = {clip(d1 / d2, 1500.0):.1f}')
    elif length_g1 < data_length or length_g2 < data_length:
        print(f'   {WHITE}UNBOUNDED{NORMAL} value = {10.0:.1f} {10.0:.1f} ratio = {-100.0:.1f} {length_g1} / {length_g2} lines out of {data_length}')
    else:
        print(f'INCORRECT DATA SIZE: found {length_g1} / {length_g2} lines, expected {data_length}')

print(f'SCAN: {argv}', file=stderr)
scan()
