#!/usr/bin/env python3
# CUT: ./model-trid.py cut 8    256   100 0 30 >/dev/null
#  NM: ./model-trid.py  nm 8 1.0e-6 10000 0 30 >/dev/null
from sys import argv
from random import random

def trid (position):
	value = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		xo_1 = position[j - 1] - j
		value += (xo - 1.0)**2 - (xo * xo_1 if j > 0 else 0.0)
	return value

if argv[1] == "cut":
	from cut import coa
	dim = int(argv[2])
	m = int(argv[3])
	max_iter = int(argv[4])
	lower = float(argv[5])
	upper = float(argv[6])
	coa(trid, dim, m, max_iter, lower, upper)
elif argv[1] == "nm":
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[2])
	tolerance = float(argv[3])
	max_iter = int(argv[4])
	lower = float(argv[5])
	upper = float(argv[6])
	nelder_mead(trid, dim, tolerance, max_iter, Simplex(trid, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
