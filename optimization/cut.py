# python implementation of cut optimization algorithm
#  (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
from sys import stderr
from math import inf, pow
from random import random

def rand_range (lower, upper):
    return (upper - lower) * random() + lower

class Point:
    def __init__(self, d):
        self.x = [0.0 for _ in range(d)]
        self.f = inf

class Population:
    def __init__(self, cost, n, m, iterations, min_x, max_x):
        self.upper = [max_x for _ in range(n)]
        self.lower = [min_x for _ in range(n)]
        self.agents = [Point(n) for _ in range(m)]
        self.evaluations = 0
        for i in range(m):
            for j in range(n):
                self.agents[i].x[j] = rand_range(min_x, max_x)
            self.agents[i].f = cost(self.agents[i].x)
            self.evaluations += 1
        self.best = self.agents[0]
        for i in range(1, m):
            if self.agents[i].f < self.best.f:
                self.best = self.agents[i]
        self.delta = pow(1.0e-6, 1.0 / iterations)

def coa(cost, n, m, max_i, min_x, max_x):
    p = Population(cost, n, m, max_i, min_x, max_x)
    iteration = 0
    while iteration < max_i:
        side = 0.5 * p.delta * (p.upper[0] - p.lower[0])
        for j in range(n):
            upper = p.best.x[j] + side
            lower = p.best.x[j] - side
            if lower < min_x:
                upper += min_x - lower
                lower = min_x
            elif upper > max_x:
                lower += max_x - upper
                upper = max_x
            p.upper[j] = upper
            p.lower[j] = lower
        for i in range(m):
            if p.agents[i] != p.best:
                for j in range(n):
                    p.agents[i].x[j] = rand_range(p.lower[j], p.upper[j])
                p.agents[i].f = cost(p.agents[i].x)
                p.evaluations += 1
        p.updated = False
        for i in range(m):
            if p.agents[i].f < p.best.f:
                p.best = p.agents[i]
                p.updated = True
        iteration += 1
        if p.updated or iteration == max_i:
            print(f'{iteration:4} {p.evaluations:6}  [ ', end='')
            print(''.join(f'{term: .6f} ' for term in p.best.x), end='')
            print(f']  {p.best.f: .6f} ')

    print(f'{iteration:4} {p.evaluations:6}  [ ', end='', file=stderr)
    print(''.join(f'{term: .6f} ' for term in p.best.x), end='', file=stderr)
    print(f']  {p.best.f: .6f} ', file=stderr)

print(__name__ + " module loaded", file=stderr)
