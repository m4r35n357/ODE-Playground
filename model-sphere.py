#!/usr/bin/env python3
from sys import argv
from random import random

def sphere (position):
	value = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		value += xo**2
	return value

if len(argv) == 6:
	from cut import coa
	dim = int(argv[1])
	m = int(argv[2])
	max_iter = int(argv[3])
	lower = float(argv[4])
	upper = float(argv[5])
	coa(sphere, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[1])
	max_iter = int(argv[2])
	lower = float(argv[3])
	upper = float(argv[4])
	nelder_mead(sphere, dim, max_iter, Simplex(sphere, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
