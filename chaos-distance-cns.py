#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv, stderr
from math import sqrt

RED = '\x1b[1;31m'
GREEN = '\x1b[1;32m'
YELLOW = '\x1b[1;33m'
WHITE = '\x1b[1;37m'
NORMAL = '\x1B[0m'

def data(line):
    d = []
    for number in line.split():
        d.append(float(number))
    return d

def separation(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)

def clip (a, limit):
    return a if a < limit else limit

def read_files(file_a, file_b):
    data_a, data_b = [], []
    with open(file_a) as a, open(file_b) as b:
        for line_a, line_b in zip(a, b):
            data_a.append(data(line_a))
            data_b.append(data(line_b))
    return data_a, data_b

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the expected data file length, and lower and uppper thresholds <<<")
    t_1, r_1 = read_files('/tmp/dataA', '/tmp/dataB')
    l_g1 = len(r_1)
    l_data, limit_threshold, chaos_threshold = int(argv[1]), float(argv[2]), float(argv[3])
    s1, max_s1 = 0.0, 0.0
    if l_g1 == l_data:
        for r, t in zip(r_1, t_1):
            s1 = separation(r, t)
            max_s1 = s1 if s1 > max_s1 else max_s1
        if max_s1 <= limit_threshold:
            print(f' {GREEN}LIMIT-CYCLE{NORMAL} value = {max_s1:.3e}')
        elif max_s1 > chaos_threshold:
            print(f'     {RED}CHAOTIC{NORMAL} value = {max_s1:.3e}')
        else:
            print(f'{YELLOW}UNCLASSIFIED{NORMAL} value = {max_s1:.3e}')
    elif l_g1 < l_data:
        print(f'   {WHITE}UNBOUNDED{NORMAL} value = {-1.0:.1f}')
    else:
        print(f'INCORRECT DATA SIZE: read {l_g1} lines, expected {l_data}')

print(f'SCAN (CNS): {argv}', file=stderr)
scan()
