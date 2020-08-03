#!/usr/bin/env python3

# Example: ./ic .00001 ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .125 .1
#          ./chaos-distance.py /tmp/dataA /tmp/dataB /tmp/dataC /tmp/dataD /tmp/dataE /tmp/dataF /tmp/dataG 50001 0.000001

from sys import argv, stderr
from math import sqrt

def line_to_data(line):
    data = []
    for number in line.split():
        data.append(float(number))
    return data

def scan():
    if len(argv) != 4:
        raise Exception(">>> ERROR! Please the intended data file length, and two separations <<<")
    data_a1, data_b1, data_c1, data_d1, data_e1, data_f1, data_g1 = [], [], [], [], [], [], []
    data_a2, data_b2, data_c2, data_d2, data_e2, data_f2, data_g2 = [], [], [], [], [], [], []
    with open('/tmp/dataA') as a, open('/tmp/dataB') as b, open('/tmp/dataC') as c, open('/tmp/dataD') as d, open('/tmp/dataE') as e, open('/tmp/dataF') as f, open('/tmp/dataG') as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a1.append(line_to_data(line_a))
            data_b1.append(line_to_data(line_b))
            data_c1.append(line_to_data(line_c))
            data_d1.append(line_to_data(line_d))
            data_e1.append(line_to_data(line_e))
            data_f1.append(line_to_data(line_f))
            data_g1.append(line_to_data(line_g))
    data_g1_length = len(data_g1)
    with open('/tmp/dataA') as a, open('/tmp/dataB') as b, open('/tmp/dataC') as c, open('/tmp/dataD') as d, open('/tmp/dataE') as e, open('/tmp/dataF') as f, open('/tmp/dataG') as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a2.append(line_to_data(line_a))
            data_b2.append(line_to_data(line_b))
            data_c2.append(line_to_data(line_c))
            data_d2.append(line_to_data(line_d))
            data_e2.append(line_to_data(line_e))
            data_f2.append(line_to_data(line_f))
            data_g2.append(line_to_data(line_g))
    data_g2_length = len(data_g2)
    assert(data_g1_length == data_g2_length)
    data_length = int(argv[1])
    separation1 = float(argv[2])
    separation2 = float(argv[3])
    results1, results2 = [], []
    r1, r2 = 0.0, 0.0
    if data_g1_length == data_length:
        for ref, trajectory in zip(data_g1, zip(data_a1, data_b1, data_c1, data_d1, data_e1, data_f1)):
            d_max = 0.0
            for point in trajectory:
                dist = sqrt((ref[0] - point[0])**2 + (ref[1] - point[1])**2 + (ref[2] - point[2])**2)
                d_max = dist if dist > d_max else d_max
            results1.append(d_max)
        r1 = results1[-1]
    if data_g2_length == data_length:
        for ref, trajectory in zip(data_g2, zip(data_a2, data_b2, data_c2, data_d2, data_e2, data_f2)):
            d_max = 0.0
            for point in trajectory:
                dist = sqrt((ref[0] - point[0])**2 + (ref[1] - point[1])**2 + (ref[2] - point[2])**2)
                d_max = dist if dist > d_max else d_max
            results2.append(d_max)
        r2 = results1[-1]
    if data_g1_length < data_length or data_g2_length < data_length:
        print(f'UNBOUNDED ({"dummy value"} = {-0.1:.1f} {data_g1_length} / {data_g2_length} lines out of {data_length})')
    elif r1 < separation1 or r2 < separation2:
        print(f'   STABLE ({"final value"} = {r1:.3e} {r2:.3e})')
    elif r1 / r2 < 0.1 * separation1 / separation2:
        print(f'  CHAOTIC ({"final value"} = {r1:.3e} {r2:.3e})')
    elif 0.9 * separation1 / separation2 < r1 / r2 < 1.1 * separation1 / separation2:
        print(f' PERIODIC ({"final value"} = {r1:.3e} {r2:.3e})')
    else:
        print(f'  UNKNOWN ')

print(f'SCAN: {argv}', file=stderr)
scan()
