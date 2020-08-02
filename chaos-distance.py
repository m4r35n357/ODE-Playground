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
    if len(argv) != 10:
        raise Exception(">>> ERROR! Please supply seven file names, the intended data file length, and a threshold <<<")
    data_a, data_b, data_c, data_d, data_e, data_f, data_g = [], [], [], [], [], [], []
    with open(argv[1]) as a, open(argv[2]) as b, open(argv[3]) as c, open(argv[4]) as d, open(argv[5]) as e, open(argv[6]) as f, open(argv[7]) as g:
        for line_a, line_b, line_c, line_d, line_e, line_f, line_g in zip(a, b, c, d, e, f, g):
            data_a.append(line_to_data(line_a))
            data_b.append(line_to_data(line_b))
            data_c.append(line_to_data(line_c))
            data_d.append(line_to_data(line_d))
            data_e.append(line_to_data(line_e))
            data_f.append(line_to_data(line_f))
            data_g.append(line_to_data(line_g))
    data_g_length = len(data_g)
    data_length = int(argv[8])
    threshold = float(argv[9])
    results = []
    r = 0.0
    if data_g_length == data_length:
        for ref, trajectory in zip(data_g, zip(data_a, data_b, data_c, data_d, data_e, data_f)):
            d_max = 0.0
            for point in trajectory:
                dist = sqrt((ref[0] - point[0])**2 + (ref[1] - point[1])**2 + (ref[2] - point[2])**2)
                d_max = dist if dist > d_max else d_max
            results.append(d_max)
        r = results[-1]
    if data_g_length < data_length:
        print(f'UNBOUNDED ({"dummy value"} = {-0.1:.3e} {data_g_length} lines out of {data_length})')
    elif r >= threshold:
        print(f'  CHAOTIC ({"final value"} = {r:.3e} >= {threshold:.3e})')
    elif r >= threshold / 1000:
        print(f' PERIODIC ({"final value"} = {r:.3e}  < {threshold:.3e})')
    else:
        print(f'   STABLE ({"final value"} = {r:.3e}  < {(threshold / 1000):.3e})')


def main():
    print(f'SCAN: {argv}', file=stderr)
    scan()

main()
