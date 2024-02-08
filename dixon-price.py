#!/usr/bin/env python3

from sys import argv
import numpy as np
from spiral_optimization import find_centre, new_point

def cost(x, dim):
    value = (x[0] - 1.0)**2
    for j in range(1, dim):
        value += (j + 1) * (2.0 * x[j]**2 - x[j - 1])**2
    return value

np.set_printoptions(precision=4, suppress=True, sign=" ")
np.set_printoptions(formatter={'float': '{: 0.4f}'.format})
np.random.seed(4)

n = int(argv[1])  # problem dimension
m = int(argv[2])  # number of points / possible solutions
max_iter = int(argv[3])
delta = float(argv[4])
r = pow(delta, 1.0 / max_iter)  # "Periodic Descent Direction Setting" rule
print("Setting number of points m = %d " % m)
print("Setting max_iter = %d " % max_iter)
print("Setting r = %0.3f " % r)

points = np.random.uniform(low=-10.0, high=10.0, size=(m, n))
centre = find_centre(points, m, cost)
print("\nInitial centre (best) point: ")
print(centre)

for itr in range(max_iter):
    if itr % 100 == 0:
        print("itr = %5d  curr centre = " % itr, end="")
        print(centre)
    for i in range(m):  # move each point towards the centre
        p = points[i]
        points[i] = new_point(p, r, centre)
    centre = find_centre(points, m, cost)  # find new centre

print("\nBest solution found: ")
print(centre, cost(centre, n))
