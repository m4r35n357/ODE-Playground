# python implementation of cut optimization algorithm
from random import random
from math import inf, pow
from sys import stderr

def rand_range (lower, upper):
	return (upper - lower) * random() + lower

class Point:
	def __init__(self, d):
		self.x = [0.0 for _ in range(d)]
		self.f = inf

class Population:
	def __init__(self, cost, n, m, iterations, min_x, max_x):
		self.max = [max_x for _ in range(n)]
		self.min = [min_x for _ in range(n)]
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
		self.updated = False

def coa(cost, n, m, max_i, min_x, max_x):
	p = Population(cost, n, m, max_i, min_x, max_x)
	iteration = 0
	while iteration < max_i:
		side = 0.5 * p.delta * (p.max[0] - p.min[0])
		for j in range(n):
			upper = p.best.x[j] + side
			lower = p.best.x[j] - side
			if lower < min_x:
				upper += min_x - lower
				lower = min_x
			elif upper > max_x:
				lower += max_x - upper
				upper = max_x
			p.max[j] = upper
			p.min[j] = lower
		for i in range(m):
			if p.agents[i] != p.best:
				for j in range(n):
					p.agents[i].x[j] = rand_range(p.min[j], p.max[j])
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
			print(''.join(f'{term: .6e} ' for term in p.best.x), end='')
			print(f']  {p.best.f: .6e} ')

print(__name__ + " module loaded", file=stderr)
