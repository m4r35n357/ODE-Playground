#!/usr/bin/env python3

from sys import argv
from math import cos, pi, exp
from cut import coa

def easom(position):
	prod = 1.0
	psum = 0.0
	for j in range(len(position)):
		xo = position[j] - (j + 1)
		prod *= cos(xo)**2
		psum -= (xo - pi)**2
	return - prod * exp(psum)

if len(argv) != 6:
	raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(easom, n, m, max_iter, lower, upper)
