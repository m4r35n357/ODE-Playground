#!/usr/bin/env python3
# CUT: ./model-rosenbrock.py 8 256 100 0 10 >/dev/null
#  NM: ./model-rosenbrock.py 8 10000 0 10 >/dev/null
from sys import argv
from random import random

def rosenbrock (position):
	value = 0.0
	for j in range(len(position) - 1):
		xo = position[j] - (j + 1)
		xo_1 = position[j + 1] - (j + 2)
		value += (1.0 - xo)**2 + 100.0 * (xo_1 - xo**2)**2
	return value

if len(argv) == 6:
	from cut import coa
	dim = int(argv[1])
	m = int(argv[2])
	max_iter = int(argv[3])
	lower = float(argv[4])
	upper = float(argv[5])
	coa(rosenbrock, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[1])
	max_iter = int(argv[2])
	lower = float(argv[3])
	upper = float(argv[4])
	nelder_mead(rosenbrock, dim, max_iter, Simplex(rosenbrock, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
