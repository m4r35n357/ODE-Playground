#!/usr/bin/env python3

from math import cos, pi, sqrt, fabs

import numpy as np

def sq_rt(x, dim):
    dim = len(x)
    value = 0.0
    for j in range(dim):
        xi = x[j]
        value += sqrt(fabs(xi))
    return value

def rosenbrock(x, dim):
    # min val = 0.0 at [1,1,1, . . 1]
    dim = len(x)
    z = 0.0
    for i in range(dim - 1):
        a = 100 * ((x[i + 1] - x[i] ** 2) ** 2)
        b = (1 - x[i]) ** 2
        z += a + b
    err = (z - 0.0) ** 2
    return err

def sphere(x, dim):
    dim = len(x)
    err = 0.0
    for i in range(dim):
        err += x[i] * x[i]
    return err

def rastrigin(x, dim):
    a = 10.0
    dim = len(x)
    value = a * dim
    for j in range(dim):
        xi = x[j]
        value += xi * xi - a * cos(2.0 * pi * xi)
    return value

def find_centre(points, m):
    n = len(points[0])  # n = dim
    best_err = sq_rt(points[0], n)
    idx = 0
    for i in range(m):
        err = sq_rt(points[i], n)
        if err < best_err:
            idx = i
            best_err = err
    return np.copy(points[idx])

def new_point(x, r, c):
    n = len(x)
    new_x = np.zeros(n)
    for k in range(n):
        new_x[k] = - (c[k] - r * (x[n - 1] - c[n - 1])) if k == 0 else c[k] - r * (x[k - 1] - c[k - 1])
    return new_x

def main():
    # 0. prepare
    np.set_printoptions(precision=4, suppress=True, sign=" ")
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format})
    np.random.seed(4)

    max_iter = 1000
    theta = np.pi / 2  # radians. small = curved, large = squared
    #    r = pow(1.0e-3, 1.0 / max_iter)  # closer to 1 = tight spiral, closer 0 = loose
    r = 0.98  # closer to 1 = tight spiral, closer 0 = loose
    m = 50  # number of points / possible solutions
    n = 3  # problem dimension

    print("\nSetting theta = %0.4f " % theta)
    print("Setting r = %0.2f " % r)
    print("Setting number of points m = %d " % m)
    print("Setting max_iter = %d " % max_iter)

    # 1. set up the Rotation matrix for n=3

    # 2. create m initial points and find initial centre (best point)
    print("\nCreating %d initial random points " % m)
    points = np.random.uniform(low=-5.0, high=5.0, size=(m, 3))

    centre = find_centre(points, m)
    print("\nInitial centre (best) point: ")
    print(centre)

    # 3. spiral points towards current centre, update centre, repeat
    print("\nUsing spiral dynamics optimization: ")
    for itr in range(max_iter):
        if itr % 100 == 0:
            print("itr = %5d  curr centre = " % itr, end="")
            print(centre)
        for i in range(m):  # move each point towards the centre
            x = points[i]
            points[i] = new_point(x, r, centre)
            #print(points)
            #input()
        centre = find_centre(points, m)  # find new centre

    # 4. show best centre found
    print("\nBest solution found: ")
    print(centre)

    print("\nEnd spiral dynamics optimization demo ")


if __name__ == "__main__":
    main()
