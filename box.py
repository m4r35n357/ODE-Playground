#!/usr/bin/env python3

from sys import argv
from math import sqrt
from whale import woa

def box(position):
	x = position[0]
	y = position[1]
	return A / (x * y * (A - 2.0 * x * y) / (2.0 * (x + y)))

model = box
if len(argv) == 1:
	dim = 2
	num_whales = 50
	max_iter = 100
else:
	dim = int(argv[1])
	num_whales = int(argv[2])
	max_iter = int(argv[3])

A = 48.0
min_edge = 0.001
best_x = woa(model, max_iter, num_whales, dim, min_edge, sqrt(0.5 * A) - min_edge)

print("Best solution found:")
for i in range(dim):
	print(f'{best_x[i]: .{6}f}', end=" ")
print("")
print("fitness of best solution = %.6f\n" % model(best_x))
