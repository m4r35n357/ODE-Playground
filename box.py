#!/usr/bin/env python3

from sys import argv
from whale import woa

def box(position):
	A = 48.0
	xi = position[0]
	yi = position[1]
	return A / (xi * yi * (A - 2.0 * xi * yi) / (2.0 * (xi + yi)))

model = box
if len(argv) == 1:
	dim = 2
	num_whales = 50
	max_iter = 100
else:
	dim = int(argv[1])
	num_whales = int(argv[2])
	max_iter = int(argv[3])

best_x = woa(model, max_iter, num_whales, dim, 0.001, 4.5)

print("Best solution found:")
print(["%.6f" % best_x[k] for k in range(dim)])
print("fitness of best solution = %.6f\n" % model(best_x))
for i in range(dim):
	print("", best_x[i], end=" ")
print("")
