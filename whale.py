# python implementation of whale optimization algorithm (WOA)
# minimizing rastrigin and sphere function
import random
import math # cos() for Rastrigin
import copy # array-copying convenience
import sys # max float

# -------fitness functions---------

# rastrigin function
def fitness_rastrigin(position):
	A = 10.0
	n = len(position)
	fitness_value = A * n
	for i in range(len(position)):
		xi = position[i]
		fitness_value += xi * xi - A * math.cos(2 * math.pi * xi)
	return fitness_value


# sphere function
def fitness_sphere(position):
	fitness_value = 0.0
	for i in range(len(position)):
		xi = position[i]
		fitness_value += (xi * xi)
	return fitness_value


# -------------------------


# whale class
class whale:
	def __init__(self, fitness, dim, minx, maxx, seed):
		self.rnd = random.Random(seed)
		self.position = [0.0 for i in range(dim)]

		for i in range(dim):
			self.position[i] = ((maxx - minx) * self.rnd.random() + minx)

		self.fitness = fitness(self.position) # curr fitness


# whale optimization algorithm(WOA)
def woa(fitness, max_iter, n, dim, minx, maxx):
	rnd = random.Random(0)

	# create n random whales
	whalePopulation = [whale(fitness, dim, minx, maxx, i) for i in range(n)]

	# compute the value of best_position and best_fitness in the whale Population
	Xbest = [0.0 for i in range(dim)]
	Fbest = sys.float_info.max

	for i in range(n): # check each whale
		if whalePopulation[i].fitness < Fbest:
			Fbest = whalePopulation[i].fitness
			Xbest = copy.copy(whalePopulation[i].position)

	# main loop of woa
	Iter = 0
	while Iter < max_iter:

		# after every 10 iterations
		# print iteration number and best fitness value so far
		if Iter % 10 == 0 and Iter > 1:
			print("Iter = " + str(Iter) + " best fitness = %.3f" % Fbest)

		# linearly decreased from 2 to 0
		a = 2 * (1 - Iter / max_iter)
		a2=-1+Iter*((-1)/max_iter)

		for i in range(n):
			A = 2 * a * rnd.random() - a
			C = 2 * rnd.random()
			b = 1
			l = (a2-1)*rnd.random()+1
			p = rnd.random()

			D = [0.0 for i in range(dim)]
			D1 = [0.0 for i in range(dim)]
			Xnew = [0.0 for i in range(dim)]
			Xrand = [0.0 for i in range(dim)]
			if p < 0.5:
				if abs(A) > 1:
					for j in range(dim):
						D[j] = abs(C * Xbest[j] - whalePopulation[i].position[j])
						Xnew[j] = Xbest[j] - A * D[j]
				else:
					p = random.randint(0, n - 1)
					while p == i:
						p = random.randint(0, n - 1)

					Xrand = whalePopulation[p].position

					for j in range(dim):
						D[j] = abs(C * Xrand[j] - whalePopulation[i].position[j])
						Xnew[j] = Xrand[j] - A * D[j]
			else:
				for j in range(dim):
					D1[j] = abs(Xbest[j] - whalePopulation[i].position[j])
					Xnew[j] = D1[j] * math.exp(b * l) * math.cos(2 * math.pi * l) + Xbest[j]

			for j in range(dim):
				whalePopulation[i].position[j] = Xnew[j]

		for i in range(n):
			# if Xnew < minx OR Xnew > maxx
			# then clip it
			for j in range(dim):
				whalePopulation[i].position[j] = max(whalePopulation[i].position[j], minx)
				whalePopulation[i].position[j] = min(whalePopulation[i].position[j], maxx)

			whalePopulation[i].fitness = fitness(whalePopulation[i].position)

			if whalePopulation[i].fitness < Fbest:
				Xbest = copy.copy(whalePopulation[i].position)
				Fbest = whalePopulation[i].fitness


		Iter += 1
	# end-while

	# returning the best solution
	return Xbest


# ----------------------------


# Driver code for rastrigin function

print("\nBegin whale optimization algorithm on rastrigin function\n")
dim = 3
fitness = fitness_rastrigin

print("Goal is to minimize Rastrigin's function in " + str(dim) + " variables")
print("Function has known min = 0.0 at (", end="")
for i in range(dim - 1):
	print("0, ", end="")
print("0)")

num_whales = 50
max_iter = 100

print("Setting num_whales = " + str(num_whales))
print("Setting max_iter = " + str(max_iter))
print("\nStarting WOA algorithm\n")

best_position = woa(fitness, max_iter, num_whales, dim, -10.0, 10.0)

print("\nWOA completed\n")
print("\nBest solution found:")
print(["%.6f" % best_position[k] for k in range(dim)])
err = fitness(best_position)
print("fitness of best solution = %.6f" % err)

print("\nEnd WOA for rastrigin\n")

print()
print()

# Driver code for Sphere function
print("\nBegin whale optimization algorithm on sphere function\n")
dim = 3
fitness = fitness_sphere

print("Goal is to minimize sphere function in " + str(dim) + " variables")
print("Function has known min = 0.0 at (", end="")
for i in range(dim - 1):
	print("0, ", end="")
print("0)")

num_whales = 50
max_iter = 100

print("Setting num_whales = " + str(num_whales))
print("Setting max_iter = " + str(max_iter))
print("\nStarting WOA algorithm\n")

best_position = woa(fitness, max_iter, num_whales, dim, -10.0, 10.0)

print("\nWOA completed\n")
print("\nBest solution found:")
print(["%.6f" % best_position[k] for k in range(dim)])
err = fitness(best_position)
print("fitness of best solution = %.6f" % err)

print("\nEnd WOA for sphere\n")
