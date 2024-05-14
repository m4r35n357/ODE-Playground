#!/usr/bin/env python3

from sys import argv
from cut import coa

def rosenbrock(position):
	value = 0.0
	for j in range(len(position) - 1):
		xo = position[j] - (j + 1)
		xo_1 = position[j + 1] - (j + 2)
		value += (1.0 - xo)**2 + 100.0 * (xo_1 - xo**2)**2
	return value

if len(argv) != 6:
	raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(rosenbrock, n, m, max_iter, lower, upper)
