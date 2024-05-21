#!/usr/bin/env python3
from sys import argv
from random import random

def trid (position):
	value = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		xo_1 = position[j - 1] - j
		value += (xo - 1.0)**2 - (xo * xo_1 if j > 0 else 0.0)
	return value

if len(argv) == 6:
	from cut import coa
	dim = int(argv[1])
	m = int(argv[2])
	max_iter = int(argv[3])
	lower = float(argv[4])
	upper = float(argv[5])
	coa(trid, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[1])
	max_iter = int(argv[2])
	lower = float(argv[3])
	upper = float(argv[4])
	nelder_mead(trid, dim, max_iter, Simplex(trid, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
