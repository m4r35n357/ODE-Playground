# spiral_optimization.py
# Anaconda3-2020.02  Python 3.7.6  NumPy 1.19.5
# Windows 10 
# minimize Rosenbrock function n=3 using spiral optimization
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


def find_center(points, m):
    # center is point with smallest error
    # note: m = len(points) (number of points / solutions)
    n = len(points[0])  # n = dim
    best_err = sq_rt(points[0], n)
    idx = 0
    for i in range(m):
        err = sq_rt(points[i], n)
        if err < best_err:
            idx = i
            best_err = err
    return np.copy(points[idx])


def move_point(x, rr, r_minus_i, center):
    # move x vector to new position, spiraling around center
    # note: matmul() automatically promotes vector to matrix
    offset = np.matmul(r_minus_i, center)  # (3,3) x (3,) = (3,)
    new_x = np.matmul(rr, x) - offset  # (3,3) x (3,) - (3,) = (3,)
    return new_x


def new_point(x, r, c):
    n = len(x)
    new_x = np.zeros(n)
    for k in range(n):
        new_x[k] = - (c[k] - r * (x[n - 1] - c[n - 1])) if k == 0 else c[k] - r * (x[k - 1] - c[k - 1])
    return new_x


def make_rab(a, b, n, theta):
    # make Rotation matrix, dim = n, a, b are 1-based
    # helper for make_R() function
    a -= 1
    b -= 1  # convert a,b to 0-based
    result = np.zeros((n, n))
    for r in range(n):
        for c in range(n):
            if r == a and c == a:
                result[r][c] = np.cos(theta)
            elif r == a and c == b:
                result[r][c] = -np.sin(theta)
            elif r == b and c == a:
                result[r][c] = np.sin(theta)
            elif r == b and c == b:
                result[r][c] = np.cos(theta)
            elif r == c:
                result[r][c] = 1.0  # unfilled diagonal elements
    return result


def make_r(n, theta):
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                result[i][j] = 1.0  # Identity
    for i in range(1, n):
        for j in range(1, i + 1):
            print("making R " + str(n - i) + " " + str(n + 1 - j))
            rr = make_rab(n - i, n + 1 - j, n, theta)
            result = np.matmul(result, rr)  # compose
    return result


def main():
    print("\nBegin spiral dynamics optimization demo ")
    print("Goal is to minimize Rosenbrock function in dim n=3")
    print("Function has known min = 0.0 at (1, 1, 1) ")

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
    print("\nSetting up hard-coded spiral Rotation matrix R ")

    ct = np.cos(theta)
    st = np.sin(theta)
    r12 = np.array([[ct, -st, 0],
                    [st, ct, 0],
                    [0, 0, 1]])

    r13 = np.array([[ct, 0, -st],
                    [0, 1, 0],
                    [st, 0, ct]])

    r23 = np.array([[1, 0, 0],
                    [0, ct, -st],
                    [0, st, ct]])

    rr = np.matmul(r23, np.matmul(r13, r12))  # compose
    # R = make_r(3, theta)  # programmatic approach

    rr = r * rr  # scale / shrink

    i3 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    r_minus_i = rr - i3

    # 2. create m initial points and
    # find curr center (best point)
    print("\nCreating %d initial random points " % m)
    points = np.random.uniform(low=-5.0, high=5.0,
                               size=(m, 3))

    center = find_center(points, m)
    print("\nInitial center (best) point: ")
    print(center)

    # 3. spiral points towards curr center,
    # update center, repeat
    print("\nUsing spiral dynamics optimization: ")
    for itr in range(max_iter):
        if itr % 100 == 0:
            print("itr = %5d  curr center = " % itr, end="")
            print(center)
        for i in range(m):  # move each pt toward center
            x = points[i]
            points[i] = move_point(x, rr, r_minus_i, center)
            #points[i] = new_point(x, r, center)
            # print(points); input()
        center = find_center(points, m)  # find new center

    # 4. show best center found
    print("\nBest solution found: ")
    print(center)

    print("\nEnd spiral dynamics optimization demo ")


if __name__ == "__main__":
    main()
