#!/usr/bin/env python3

from sys import argv
from math import fsum, sqrt
from cut import coa

if len(argv) != 6:
	raise Exception('>>> Wrong number of arguments <<<')
n = int(argv[1])
m = int(argv[2])
max_iter = int(argv[3])
lower = float(argv[4])
upper = float(argv[5])

coa(lambda x: fsum(sqrt(abs(x[j] - (j + 1))) for j in range(len(x))), n, m, max_iter, lower, upper)