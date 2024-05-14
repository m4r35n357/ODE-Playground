#!/usr/bin/env python3

from sys import argv
from cut import coa

def dixon_price(position):
    value = (position[0] - 1.0)**2
    for j in range(1, len(position)):
        value += (j + 1) * (2.0 * position[j]**2 - position[j - 1])**2
    return value

if len(argv) != 6:
    raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(dixon_price, n, m, max_iter, lower, upper)
