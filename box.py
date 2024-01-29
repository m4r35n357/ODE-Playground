#!/usr/bin/env python3

from sys import argv
from math import sqrt
from whale import woa

def z(x, y):
	return (A - 2.0 * x * y) / (2.0 * (x + y))

def box(position):
	x = position[0]
	y = position[1]
	return A / (x * y * z(x, y))

model = box
dim = 2
if len(argv) == 1:
	num_whales = 50
	max_iter = 100
else:
	num_whales = int(argv[1])
	max_iter = int(argv[2])

A = 48.0
min_edge = 0.001
best_x = woa(model, max_iter, num_whales, dim, min_edge, sqrt(0.5 * A) - min_edge)

print("Best solution found:")
print(f'x {best_x[0]: .{6}f}  y {best_x[1]: .{6}f}  z {z(best_x[0], best_x[1]): .{6}f} ', end=" ")
print(f'V {best_x[0] * best_x[1] * z(best_x[0], best_x[1]): .{6}f}')
print("fitness of best solution = %.6f\n" % model(best_x))
