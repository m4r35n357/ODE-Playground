#!/usr/bin/env python3

from sys import argv
from math import cos, pi

from whale import woa

def rastrigin(position):
	a = 10.0
	n = len(position)
	value = a * n
	for j in range(n):
		xi = position[j]
		value += xi * xi - a * cos(2.0 * pi * xi)
	return value

model = rastrigin
if len(argv) == 1:
	dim = 3
	num_whales = 50
	max_iter = 100
else:
	dim = int(argv[1])
	num_whales = int(argv[2])
	max_iter = int(argv[3])

best_x = woa(model, max_iter, num_whales, dim, -10.0, 10.0)

print("Best solution found:")
print(["%.6f" % best_x[k] for k in range(dim)])
print("fitness of best solution = %.6f\n" % model(best_x))
for i in range(dim):
	print("", best_x[i], end=" ")
print("")

