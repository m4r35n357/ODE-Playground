#!/usr/bin/env python3

from sys import stderr
from math import pi, sqrt

from playground import analyze

r2 = 10.0
q = r2 / 2.0
rtq = sqrt(q - 1.0)


def playground(a, value):
    # Example: ./grtwins.py 2 1.1 1.9 1001 0 1e-9 1e-9 | ./plotMany.py 2 10
    return pi * (1.0 + a.cos) ** 1.5 * q ** 1.5 - (
            q * rtq * (a + a.sin) + 2.0 * a * rtq + 2.0 * ((rtq + (a / 2.0).tan) / (rtq - (a / 2.0).tan)).ln) - value


for result in analyze(playground, max_it=100):
    print(result, file=stderr)
