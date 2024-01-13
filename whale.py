# python implementation of whale optimization algorithm (WOA)
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
	whale = [Whale(cost, n, min_x, max_x, j) for j in range(n_whales)]
	# compute the value of best_position and best_fitness in the whale Population
	x_best = [0.0 for _ in range(n)]
	f_best = float_info.max

	for k in range(n_whales): # check each whale
		if whale[k].value < f_best:
			f_best = whale[k].value
			x_best = copy(whale[k].position)
	# main loop of woa
	iteration = 0
	while iteration < i_max:
		# after every 10 iterations print iteration number and best fitness value so far
		if iteration % 10 == 0 and iteration > 1:
			print("Iter = " + str(iteration) + " best fitness = %.3f" % f_best)
		# linearly decreased from 2 to 0
		a = 2 * (1 - iteration / i_max)
		a2 = -1 + iteration * ((-1) / i_max)
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
						d[j] = abs(cc * x_best[j] - whale[k].position[j])
						x_new[j] = x_best[j] - aa * d[j]
				else:
					p = randint(0, n_whales - 1)
					while p == k:
						p = randint(0, n_whales - 1)

					x_rand = whale[p].position

					for j in range(n):
						d[j] = abs(cc * x_rand[j] - whale[k].position[j])
						x_new[j] = x_rand[j] - aa * d[j]
			else:
				for j in range(n):
					d_1[j] = abs(x_best[j] - whale[k].position[j])
					x_new[j] = d_1[j] * exp(b * l) * cos(2 * pi * l) + x_best[j]
			for j in range(n):
				whale[k].position[j] = x_new[j]
		for k in range(n_whales):
			# if Xnew < minx OR Xnew > maxx then clip it
			for j in range(n):
				whale[k].position[j] = max(whale[k].position[j], min_x)
				whale[k].position[j] = min(whale[k].position[j], max_x)
			whale[k].value = cost(whale[k].position)
			if whale[k].value < f_best:
				x_best = copy(whale[k].position)
				f_best = whale[k].value
		iteration += 1
	# end-while
	return x_best
