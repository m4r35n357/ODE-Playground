#!/usr/bin/env python3

from spiral_optimization import soa

def cost(x, dim):
    value = (x[0] - 1.0)**2
    for j in range(1, dim):
        value += (j + 1) * (2.0 * x[j]**2 - x[j - 1])**2
    return value

soa()