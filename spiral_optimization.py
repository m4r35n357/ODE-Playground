from sys import argv, stderr
from math import cos, pi, sqrt, fabs
import numpy as np

def sq_rt(x, dim):
    dim = len(x)
    value = 0.0
    for j in range(dim):
        xi = x[j]
        value += sqrt(fabs(xi))
    return value

def dixon_price(x, dim):
    dim = len(x)
    value = (x[0] - 1.0)**2
    for j in range(1, dim):
        value += (j + 1) * (2.0 * x[j]**2 - x[j - 1])**2
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

def find_centre(points, m, function):
    n = len(points[0])  # n = dim
    best_err = function(points[0], n)
    idx = 0
    for i in range(m):
        err = function(points[i], n)
        if err < best_err:
            idx = i
            best_err = err
    centre = points[idx]
    return np.copy(centre)

def new_point(x, r, c):
    n = len(x)
    new_x = np.zeros(n)
    for k in range(n):
        new_x[k] = c[k] - r * (x[n - 1] - c[n - 1]) if k == 0 else c[k] + r * (x[k - 1] - c[k - 1])
    return new_x

def main():
    function = dixon_price

    # 0. prepare
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

    # 1. set up the Rotation matrix for n=3 - N/A - using newer method

    # 2. create m initial points and find initial centre (best point)
    points = np.random.uniform(low=-10.0, high=10.0, size=(m, n))

    centre = find_centre(points, m, function)
    print("\nInitial centre (best) point: ")
    print(centre)

    # 3. spiral points towards current centre, update centre, repeat
    for itr in range(max_iter):
        if itr % 100 == 0:
            print("itr = %5d  curr centre = " % itr, end="")
            print(centre)
        for i in range(m):  # move each point towards the centre
            x = points[i]
            points[i] = new_point(x, r, centre)
        centre = find_centre(points, m, function)  # find new centre

    # 4. show best centre found
    print("\nBest solution found: ")
    print(centre, function(centre, n))


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
