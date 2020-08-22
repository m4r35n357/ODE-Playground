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

def separation(reference, trajectories):
    s = 0.0
    for t in trajectories:
        d = sqrt((reference[0] - t[0])**2 + (reference[1] - t[1])**2 + (reference[2] - t[2])**2)
        s = d if d > s else s
    return s

def clip (a, limit):
    return a if a < limit else limit

def read_files(file_a, file_b, file_c, file_d, file_e, file_f, file_g):
    data_a, data_b, data_c, data_d, data_e, data_f, data_g = [], [], [], [], [], [], []
    with open(file_a) as a, open(file_b) as b, open(file_c) as c, open(file_d) as d, open(file_e) as e, open(file_f) as f, open(file_g) as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a.append(data(line_a))
            data_b.append(data(line_b))
            data_c.append(data(line_c))
            data_d.append(data(line_d))
            data_e.append(data(line_e))
            data_f.append(data(line_f))
            data_g.append(data(line_g))
    return (data_a, data_b, data_c, data_d, data_e, data_f), data_g

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the expected data file length, and two initial separation values (greater first!) <<<")
    t_1, r_1 = read_files('/tmp/dataA1', '/tmp/dataB1', '/tmp/dataC1', '/tmp/dataD1', '/tmp/dataE1', '/tmp/dataF1', '/tmp/dataG1')
    t_2, r_2 = read_files('/tmp/dataA2', '/tmp/dataB2', '/tmp/dataC2', '/tmp/dataD2', '/tmp/dataE2', '/tmp/dataF2', '/tmp/dataG2')
    l_g1, l_g2 = len(r_1), len(r_2)
    l_data, δ_1, δ_2 = int(argv[1]), float(argv[2]), float(argv[3])
    ratio = δ_1 / δ_2
    s1, s2, max_s1, max_s2 = 0.0, 0.0, 0.0, 0.0
    if l_g1 == l_data and l_g2 == l_data:
        for r, t in zip(r_1, t_1):
            s1 = separation(r, t)
            max_s1 = s1 if s1 > max_s1 else max_s1
        for r, t in zip(r_2, t_2):
            s2 = separation(r, t)
            max_s2 = s2 if s2 > max_s2 else max_s2
        if s1 < δ_1 and s2 < δ_2:
            print(f'   {BLUE}CONVERGED{NORMAL} value = {s1:.3e} {s2:.3e} ratio = {2.0*ratio:.1f}')
        elif 0.8 * ratio < max_s1 / max_s2 < 1.2 * ratio:
            print(f' {GREEN}LIMIT-CYCLE{NORMAL} value = {max_s1:.3e} {max_s2:.3e} ratio = {max_s1/max_s2:.1f}')
        elif 0.5 < max_s1 / max_s2 < 2.0:
            print(f'     {RED}CHAOTIC{NORMAL} value = {max_s1:.3e} {max_s2:.3e} ratio = {max_s1/max_s2:.1f}')
        else:
            print(f'{YELLOW}UNCLASSIFIED{NORMAL} value = {max_s1:.3e} {max_s2:.3e} ratio = {clip(max_s1/max_s2, 1.5*ratio):.1f}')
    elif l_g1 < l_data or l_g2 < l_data:
        print(f'   {WHITE}UNBOUNDED{NORMAL} value = {ratio:.1f} {ratio:.1f} ratio = {-0.5*ratio:.1f} {l_g1}/{l_g2} of {l_data}')
    else:
        print(f'INCORRECT DATA SIZE: read {l_g1}/{l_g2} lines, expected {l_data}')

print(f'SCAN: {argv}', file=stderr)
scan()
