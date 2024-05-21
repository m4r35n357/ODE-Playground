#!/usr/bin/env python3
from sys import argv
from math import cos, pi
from random import random

def rastrigin (position):
	a = 10.0
	n = len(position)
	value = a * n
	for j in range(n):
		xo = position[j] - (j + 1)
		value += xo * xo - a * cos(2.0 * pi * xo)
	return value

if len(argv) == 6:
	from cut import coa
	dim = int(argv[1])
	m = int(argv[2])
	max_iter = int(argv[3])
	lower = float(argv[4])
	upper = float(argv[5])
	coa(rastrigin, dim, m, max_iter, lower, upper)
elif len(argv) == 5:
	from nelder_mead import nelder_mead, Simplex
	dim = int(argv[1])
	max_iter = int(argv[2])
	lower = float(argv[3])
	upper = float(argv[4])
	nelder_mead(rastrigin, dim, max_iter, Simplex(rastrigin, [((upper - lower) * random() + lower) for _ in range(dim)], 1.0))
else:
	raise Exception('>>> Wrong number of arguments <<<')
