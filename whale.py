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
	whales = [Whale(cost, n, min_x, max_x, j) for j in range(n_whales)]
	# compute the value of best_position and best_fitness in the whale Population
	x_best = [0.0 for _ in range(n)]
	f_best = float_info.max
	for k in range(n_whales): # find best whale
		if whales[k].value < f_best:
			f_best = whales[k].value
			x_best = copy(whales[k].position)
	# main loop of woa
	iteration = 0
	while iteration < i_max:
		# after every 10 iterations print iteration number and best fitness value so far
		if iteration % 10 == 0 and iteration > 1:
			print("Iter = " + str(iteration) + " best fitness = %.3f" % f_best)
		# linearly decreased from 2 to 0
		a = 2.0 * (1 - iteration / i_max)
		a2 = -1.0 + iteration * ((-1) / i_max)
		for k in range(n_whales):
			aa = 2 * a * rnd.random() - a
			cc = 2 * rnd.random()
			b = 1
			l = (a2 - 1) * rnd.random() + 1
			p = rnd.random()
			d = [0.0 for _ in range(n)]
			d_1 = [0.0 for _ in range(n)]
			x_new = [0.0 for _ in range(n)]
			if p < 0.5:
				if abs(aa) > 1:
					for j in range(n):
						d[j] = abs(cc * x_best[j] - whales[k].position[j])
						x_new[j] = x_best[j] - aa * d[j]
				else:
					p = randint(0, n_whales - 1)
					while p == k:
						p = randint(0, n_whales - 1)

					x_rand = whales[p].position

					for j in range(n):
						d[j] = abs(cc * x_rand[j] - whales[k].position[j])
						x_new[j] = x_rand[j] - aa * d[j]
			else:
				for j in range(n):
					d_1[j] = abs(x_best[j] - whales[k].position[j])
					x_new[j] = d_1[j] * exp(b * l) * cos(2 * pi * l) + x_best[j]
			for j in range(n):
				whales[k].position[j] = x_new[j]
		for k in range(n_whales):
			# if Xnew < minx OR Xnew > maxx then clip it
			for j in range(n):
				whales[k].position[j] = max(whales[k].position[j], min_x)
				whales[k].position[j] = min(whales[k].position[j], max_x)
			whales[k].value = cost(whales[k].position)
			if whales[k].value < f_best:
				x_best = copy(whales[k].position)
				f_best = whales[k].value
		iteration += 1
	# end-while
	return x_best
