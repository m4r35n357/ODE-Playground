#!/usr/bin/env python3

"""
Example: ./grtwins.py 2 -2 2 1001 0 1e-9 1e-9 | ./plotMany.py 2 300 >/dev/null

see https://www.mathpages.com/rr/s6-05/6-05.htm, equations (7) and (8)
"""

from sys import stderr
from math import pi, sqrt, sin, cos
from playground import Solver, analyze

n = 1
r2 = 10.0
q = r2 / 2.0
rtq = sqrt(q - 1.0)


def playground(a, value):
    return pi * n * (1.0 + a.cos) ** 1.5 * q ** 1.5 - (
            q * rtq * (a + a.sin) + 2.0 * a * rtq + 2.0 * ((rtq + (a / 2.0).tan) / (rtq - (a / 2.0).tan)).ln) - value


def ratio(a):
    r1 = (1 + cos(a)) * q
    return r1, ((a + sin(a)) * sqrt(q)) / (pi * n * (1 + cos(a)) * sqrt(r1 - 3.0))


for result in analyze(playground, max_it=100):
    if result.mode == Solver.ROOT.name:
        print(result, file=stderr)
        answer = ratio(result.x)
        print(answer, file=stderr)
        print("n = {}".format(n), file=stderr)
        print("alpha = {}".format(result.x), file=stderr)
        print("r1 = {}".format(answer[0]), file=stderr)
        print("r2 = {}".format(r2), file=stderr)
        print("ratio = {}".format(answer[1]), file=stderr)
