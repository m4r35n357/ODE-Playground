# python implementation of whale optimization algorithm (WOA)
# https://www.geeksforgeeks.org/implementation-of-whale-optimization-algorithm/
from random import Random, randint
from math import exp, cos, pi
from copy import copy # array-copying convenience
from sys import float_info # max float

class Whale:
	def __init__(self, cost, n_dim, min_x, max_x, seed):
		self.rnd = Random(seed)
		self.x = [0.0 for _ in range(n_dim)]
		for j in range(n_dim):
			self.x[j] = ((max_x - min_x) * self.rnd.random() + min_x)
		self.value = cost(self.x) # curr fitness

def woa(cost, i_max, n_whales, dim, min_x, max_x):
	rnd = Random(0)
	whales = [Whale(cost, dim, min_x, max_x, i) for i in range(n_whales)]
	x_p = [0.0 for _ in range(dim)]
	f_best = float_info.max
	for i in range(n_whales): # find best whale
		if whales[i].value < f_best:
			f_best = whales[i].value
			x_p = copy(whales[i].x)
	iteration = 0
	while iteration < i_max:
		if iteration % 10 == 0 and iteration > 1:
			print("Iter = " + str(iteration) + " best fitness = %.3f" % f_best)
		a = 2.0 * (1.0 - iteration / i_max)
		for i in range(n_whales):
			aa = a * (2.0 * rnd.random() - 1.0)
			cc = 2.0 * rnd.random()
			b = 1.0
			l = 2.0 * rnd.random() - 1.0
			x_next = [0.0 for _ in range(dim)]
			if rnd.random() < 0.5:
				if abs(aa) < 1.0:  # "encircling" update (1)
					for j in range(dim):
						x_next[j] = x_p[j] - aa * abs(cc * x_p[j] - whales[i].x[j])
				else:  # "searching/random" update (9)
					p = randint(0, n_whales - 1)
					while p == i:
						p = randint(0, n_whales - 1)
					x_rand = whales[p].x
					for j in range(dim):
						x_next[j] = x_rand[j] - aa * abs(cc * x_rand[j] - whales[i].x[j])
			else:  # "spiral" update (7)
				for j in range(dim):
					x_next[j] = abs(x_p[j] - whales[i].x[j]) * exp(b * l) * cos(2.0 * pi * l) + x_p[j]
			for j in range(dim):
				whales[i].x[j] = x_next[j]
		for i in range(n_whales):
			for j in range(dim):
				whales[i].x[j] = max(whales[i].x[j], min_x)
				whales[i].x[j] = min(whales[i].x[j], max_x)
			whales[i].value = cost(whales[i].x)
			if whales[i].value < f_best:
				x_p = copy(whales[i].x)
				f_best = whales[i].value
		iteration += 1
	return x_p
