#!/usr/bin/env python3

# Example: ./ic .00001 ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .125 .1
#          ./chaos-distance.py /tmp/dataA /tmp/dataB /tmp/dataC /tmp/dataD /tmp/dataE /tmp/dataF /tmp/dataG 50001 0.000001

from sys import argv, stderr
from math import sqrt, log

def line_to_data(line):
    data = []
    for number in line.split():
        data.append(float(number))
    return data

def distance(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the intended data file length, and two separations <<<")
    # read first set of data
    data_a1, data_b1, data_c1, data_d1, data_e1, data_f1, data_g1 = [], [], [], [], [], [], []
    with open('/tmp/dataA1') as a, open('/tmp/dataB1') as b, open('/tmp/dataC1') as c, open('/tmp/dataD1') as d, open('/tmp/dataE1') as e, open('/tmp/dataF1') as f, open('/tmp/dataG1') as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a1.append(line_to_data(line_a))
            data_b1.append(line_to_data(line_b))
            data_c1.append(line_to_data(line_c))
            data_d1.append(line_to_data(line_d))
            data_e1.append(line_to_data(line_e))
            data_f1.append(line_to_data(line_f))
            data_g1.append(line_to_data(line_g))
    data_g1_length = len(data_g1)
    # read second set of data
    data_a2, data_b2, data_c2, data_d2, data_e2, data_f2, data_g2 = [], [], [], [], [], [], []
    with open('/tmp/dataA2') as a, open('/tmp/dataB2') as b, open('/tmp/dataC2') as c, open('/tmp/dataD2') as d, open('/tmp/dataE2') as e, open('/tmp/dataF2') as f, open('/tmp/dataG2') as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a2.append(line_to_data(line_a))
            data_b2.append(line_to_data(line_b))
            data_c2.append(line_to_data(line_c))
            data_d2.append(line_to_data(line_d))
            data_e2.append(line_to_data(line_e))
            data_f2.append(line_to_data(line_f))
            data_g2.append(line_to_data(line_g))
    data_g2_length = len(data_g2)
    # process data sets
    data_length = int(argv[1])
    separation1, separation2 = float(argv[2]), float(argv[3])
    slope = separation1 / separation2
    r1, r2 = 0.0, 0.0
    if data_g1_length == data_length and data_g2_length == data_length:
        for ref, trajectory in zip(data_g1, zip(data_a1, data_b1, data_c1, data_d1, data_e1, data_f1)):
            d_max = 0.0
            for point in trajectory:
                dist = distance(ref, point)
                d_max = dist if dist > d_max else d_max
            r1 = d_max
        for ref, trajectory in zip(data_g2, zip(data_a2, data_b2, data_c2, data_d2, data_e2, data_f2)):
            d_max = 0.0
            for point in trajectory:
                dist = distance(ref, point)
                d_max = dist if dist > d_max else d_max
            r2 = d_max
    # analyze data
    if data_g1_length < data_length or data_g2_length < data_length:
        print(f'  UNBOUNDED ({"dummy values"} = {-0.1:.1f} {data_g1_length} / {data_g2_length} lines out of {data_length})')
    elif r1 < separation1 and r2 < separation2:
        print(f'  CONVERGED ({"final values"} = {r1:.3e} < {separation1:.1e}, {r2:.3e} < {separation2:.1e})')
    elif 0.9 * slope < r1 / r2 < 1.1 * slope:
        print(f'LIMIT CYCLE ({"final values"} = {r1:.3e} {r2:.3e} ratio = {r1 / r2:.1e})')
    elif 0.5 < r1 / r2 < 2.0:
        print(f'    CHAOTIC ({"final values"} = {r1:.3e} {r2:.3e} ratio = {r1 / r2:.1e})')
    else:
        print(f'    UNKNOWN ({"final values"} = {r1:.3e} {r2:.3e} ratio = {r1 / r2:.1e} {10.0 * log(r1 / r2):.1e})')

print(f'SCAN: {argv}', file=stderr)
scan()
