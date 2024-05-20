from math import sqrt
from sys import stderr, argv
from copy import copy
from operator import attrgetter

class Point:
    def __init__(self, n):
        self.x = [0.0 for _ in range(n)]
        self.f = 0.0

class Simplex:
    def __init__(self, function, x0, size):
        n = len(x0)
        self.iterations = self.evaluations = 0
        self.vertices = [Point(n) for _ in range(n + 1)]
        b = 0.0
        for j in range(n):
            c = sqrt(1.0 - b)
            self.vertices[j].x[j] = c
            r = - (1.0 / n + b) / c
            for i in range(j + 1, n + 1):
                self.vertices[i].x[j] = r
            b += r * r
        for i in range(n + 1):
            for j in range(n):
                self.vertices[i].x[j] = size * self.vertices[i].x[j] + x0[j]
            self.vertices[i].f = function(self.vertices[i].x)
            self.evaluations += 1
        self.vertices.sort(key=attrgetter('f'))

def project(start, factor, pa, pb):
    return [start.x[j] + factor * (pb.x[j] - pa.x[j]) for j in range(len(start))]

def nelder_mead(function, n, max_iterations, s):
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5
    centroid = [0.0 for _ in range(n)]
    for j in range(n):
        centroid[j] = 0.0
        for i in range(n):
            centroid[j] += s.vertices[i].x[j]
    reflect = Point(n)
    reflect.x = [centroid[j] + alpha * (s.vertices[n].x[j] - centroid[j]) for j in range(n)]
    reflect.f = function(reflect.x)
    s.evaluations += 1
    while s.iterations < max_iterations:
        shrink = False
        if s.vertices[0].f <= reflect.f < s.vertices[n - 1].f:
            s.vertices[n] = copy(reflect)
        elif reflect.f < s.vertices[0].f:
            expand = Point(n)
            expand.x = [centroid[j] + gamma * (s.vertices[n].x[j] - centroid[j]) for j in range(n)]
            expand.f = function(expand.x)
            s.evaluations += 1
            if expand.f < reflect.f:
                s.vertices[n] = copy(expand)
            else:
                s.vertices[n] = copy(reflect)
        elif reflect.f < s.vertices[n].f:
            contract_out = Point(n)
            contract_out.x = [centroid[j] + rho * (s.vertices[n].x[j] - centroid[j]) for j in range(n)]
            contract_out.f = function(contract_out.x)
            s.evaluations += 1
            if contract_out.f < reflect.f:
                shrink = True
        else:
            contract_in = Point(n)
            contract_in.x = [centroid[j] - rho * (s.vertices[n].x[j] - centroid[j]) for j in range(n)]
            contract_in.f = function(contract_in.x)
            s.evaluations += 1
            if contract_in.f < reflect.f:
                shrink = True
        if shrink:
            for i in range(n + 1):
                for j in range(n):
                    s.vertices[i].x[j] = [s.vertices[0].x[j] - sigma * (s.vertices[i].x[j] - s.vertices[0].x[j]) for j in range(n)]
                s.vertices[i].f = function(s.vertices[i].x)
                s.evaluations += 1
        s.vertices.sort(key=attrgetter('f'))
        s.iterations += 1
        print(f'{s.iterations:4} {s.evaluations:6}  [ ', end='')
        print(''.join(f'{term: .6f} ' for term in s.vertices[0].x), end='')
        print(f']  {s.vertices[0].f: .6f} ')

def rosenbrock(position):
    value = 0.0
    for j in range(len(position) - 1):
        xi = position[j]
        value += (1.0 - xi)**2 + 100.0 * (position[j + 1] - xi * xi)**2
    return value

if __name__ == "__main__":
    dim = int(argv[1])
    max_iter = int(argv[2])
    initial_point = [0.0, 0.0, 0.0]
    initial_simplex = Simplex(rosenbrock, initial_point, 1.0)
    nelder_mead(rosenbrock, dim, max_iter, initial_simplex)
else:
    print(__name__ + "nelder_mead module loaded", file=stderr)
