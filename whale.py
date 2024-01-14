# python implementation of whale optimization algorithm (WOA)
# https://www.geeksforgeeks.org/implementation-of-whale-optimization-algorithm/
from random import Random, randint
from math import exp, cos, pi
from copy import copy # array-copying convenience
from sys import float_info # max float

class Whale:
	def __init__(self, cost, n, min_x, max_x, seed):
		self.rnd = Random(seed)
		self.position = [0.0 for _ in range(n)]
		for j in range(n):
			self.position[j] = ((max_x - min_x) * self.rnd.random() + min_x)
		self.value = cost(self.position) # curr fitness

def woa(cost, i_max, n_whales, n, min_x, max_x):
	rnd = Random(0)
	# create n random whales
	whales = [Whale(cost, n, min_x, max_x, i) for i in range(n_whales)]
	# compute the value of best_position and best_fitness in the whale Population
	x_p = [0.0 for _ in range(n)]
	f_best = float_info.max
	for i in range(n_whales): # find best whale
		if whales[i].value < f_best:
			f_best = whales[i].value
			x_p = copy(whales[i].position)
	# main loop of woa
	iteration = 0
	while iteration < i_max:
		if iteration % 10 == 0 and iteration > 1:
			print("Iter = " + str(iteration) + " best fitness = %.3f" % f_best)
		a = 2.0 * (1.0 - iteration / i_max)
		for i in range(n_whales):
			aa = a * (2.0 * rnd.random() - 1.0)
			cc = 2.0 * rnd.random()
			b = 1.0
			l = 1.0 - (iteration / i_max + 2.0) * rnd.random()
#			l = 2.0 * rnd.random() - 1.0
			p = rnd.random()
			d = [0.0 for _ in range(n)]
			d_dash = [0.0 for _ in range(n)]
			d_ddash = [0.0 for _ in range(n)]
			x_new = [0.0 for _ in range(n)]
			if p < 0.5:
				if abs(aa) > 1.0:
					for j in range(n):
						d[j] = abs(cc * x_p[j] - whales[i].position[j])
						x_new[j] = x_p[j] - aa * d[j]
				else:
					p = randint(0, n_whales - 1)
					while p == i:
						p = randint(0, n_whales - 1)
					x_rand = whales[p].position
					for j in range(n):
						d_ddash[j] = abs(cc * x_rand[j] - whales[i].position[j])
						x_new[j] = x_rand[j] - aa * d_ddash[j]
			else:
				for j in range(n):
					d_dash[j] = abs(x_p[j] - whales[i].position[j])
					x_new[j] = d_dash[j] * exp(b * l) * cos(2 * pi * l) + x_p[j]
			for j in range(n):
				whales[i].position[j] = x_new[j]
		for i in range(n_whales):
			for j in range(n):
				whales[i].position[j] = max(whales[i].position[j], min_x)
				whales[i].position[j] = min(whales[i].position[j], max_x)
			whales[i].value = cost(whales[i].position)
			if whales[i].value < f_best:
				x_p = copy(whales[i].position)
				f_best = whales[i].value
		iteration += 1
	return x_p
