#!/usr/bin/env python3

from sys import argv
from cut import coa

def trid(position):
	value = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		xo_1 = position[j - 1] - j
		value += (xo - 1.0)**2 - (xo * xo_1 if j > 0 else 0.0)
	return value

if len(argv) != 6:
	raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(trid, n, m, max_iter, lower, upper)
