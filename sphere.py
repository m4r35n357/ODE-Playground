#!/usr/bin/env python3

from whale import woa

def sphere(position):
	value = 0.0
	for j in range(len(position)):
		xi = position[j]
		value += xi * xi
	return value

print("\nBegin whale optimization algorithm on sphere function\n")
dim = 3
model = sphere

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

best_x = woa(model, max_iter, num_whales, dim, -10.0, 10.0)

print("\nWOA completed")
print("Best solution found:")
print(["%.6f" % best_x[k] for k in range(dim)])
print("fitness of best solution = %.6f\n" % model(best_x))
for i in range(dim):
	print("", best_x[i], end=" ")
print("")
