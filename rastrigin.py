#!/usr/bin/env python3

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

print("\nBegin whale optimization algorithm on rastrigin function\n")
dim = 3
model = rastrigin

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

best_x = woa(model, max_iter, num_whales, dim, -10.0, 10.0)

print("\nWOA completed")
print("Best solution found:")
print(["%.6f" % best_x[k] for k in range(dim)])
err = model(best_x)
print("fitness of best solution = %.6f\n" % err)
