#!/usr/bin/env python3
# CUT: ./model-easom.py cut 8    256   100 0 10 >/dev/null
#  NM: ./model-easom.py  nm 8 1.0e-6 10000 0 10 >/dev/null
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

if argv[1] == "cut":
	from cut import coa
	dim = int(argv[2])
	m = int(argv[3])
	max_iter = int(argv[4])
	lower = float(argv[5])
	upper = float(argv[6])
	coa(easom, dim, m, max_iter, lower, upper)
elif argv[1] == "nm":
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[2])
	tolerance = float(argv[3])
	max_iter = int(argv[4])
	lower = float(argv[5])
	upper = float(argv[6])
	nelder_mead(easom, dim, tolerance, max_iter, Simplex(easom, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
