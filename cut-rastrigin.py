#!/usr/bin/env python3

from sys import argv
from math import cos, pi
from cut import coa

def rastrigin(position):
	a = 10.0
	n = len(position)
	value = a * n
	for j in range(n):
		xo = position[j] - (j + 1)
		value += xo * xo - a * cos(2.0 * pi * xo)
	return value

if len(argv) != 6:
	raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(rastrigin, n, m, max_iter, lower, upper)
