#!/usr/bin/env python3
# CUT: ./model-easom.py 8 256 100 0 10 >/dev/null
#  NM: ./model-easom.py 8 10000 0 10 >/dev/null
from sys import argv
from math import cos, pi, exp
from random import random

def easom (position):
	prod = 1.0
	psum = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		prod *= cos(xo)**2
		psum -= (xo - pi)**2
	return - prod * exp(psum)

if len(argv) == 6:
	from cut import coa
	dim = int(argv[1])
	m = int(argv[2])
	max_iter = int(argv[3])
	lower = float(argv[4])
	upper = float(argv[5])
	coa(easom, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[1])
	max_iter = int(argv[2])
	lower = float(argv[3])
	upper = float(argv[4])
	nelder_mead(easom, dim, max_iter, Simplex(easom, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
