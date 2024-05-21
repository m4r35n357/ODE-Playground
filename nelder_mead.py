# Nelder-Mead Algorithm in Python
from sys import stderr
from math import sqrt, fsum
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
        self.centroid = Point(n)
        self.reflect = Point(n)
        self.expand = Point(n)
        self.contract_out = Point(n)
        self.contract_in = Point(n)
        b = 0.0
        for j in range(n):
            c = sqrt(1.0 - b)
            self.vertices[j].x[j] = c
            r = - (1.0 / n + b) / c
            for i in range(j + 1, n + 1):
                self.vertices[i].x[j] = r
            b += r * r
        for i in range(n + 1):
            self.vertices[i].x = [size * self.vertices[i].x[j] + x0[j] for j in range(n)]
            self.vertices[i].f = function(self.vertices[i].x)
            self.evaluations += 1
        sort(self, n)

def project(new, s, function, n, factor, pa, pb):
    new.x = [pb.x[j] + factor * (pb.x[j] - pa.x[j]) for j in range(n)]
    new.f = function(new.x)
    s.evaluations += 1
    return new

def distance (a, b, n):
    return sqrt(fsum((a.x[i] - b.x[i])**2 for i in range(n)))

def sort (s, n):
    s.vertices.sort(key=attrgetter('f'))
    s.centroid.x = [fsum(s.vertices[i].x[j] for i in range(n)) / n for j in range(n)]

def nelder_mead(function, n, max_iterations, s):
    while ((distance(s.vertices[0], s.vertices[n], n) > 1.0e-6) or (s.vertices[0].f - s.vertices[n].f > 1.0e-6)) and s.iterations < max_iterations:
        shrink = False
        s.reflect = project(s.reflect, s, function, n, 1.0, s.vertices[n], s.centroid)
        if s.vertices[0].f <= s.reflect.f < s.vertices[n - 1].f:
            s.vertices[n] = copy(s.reflect)
        elif s.reflect.f < s.vertices[0].f:
            s.expand = project(s.expand, s, function, n, 2.0, s.vertices[n], s.centroid)
            if s.expand.f < s.reflect.f:
                s.vertices[n] = copy(s.expand)
            else:
                s.vertices[n] = copy(s.reflect)
        elif s.reflect.f < s.vertices[n].f:
            s.contract_out = project(s.contract_out, s, function, n, 0.5, s.vertices[n], s.centroid)
            if s.contract_out.f < s.reflect.f:
                s.vertices[n] = copy(s.contract_out)
            else:
                shrink = True
        else:
            s.contract_in = project(s.contract_in, s, function, n, -0.5, s.vertices[n], s.centroid)
            if s.contract_in.f < s.vertices[n].f:
                s.vertices[n] = copy(s.contract_in)
            else:
                shrink = True
        if shrink:
            for i in range(1, n + 1):
                s.vertices[i] = project(s.vertices[i], s, function, n, -0.5, s.vertices[i], s.vertices[0])
        sort(s, n)
        s.iterations += 1
        print(f'{s.iterations:4} {s.evaluations:6}  [ ', end='')
        print(''.join(f'{term: .6f} ' for term in s.vertices[0].x), end='')
        print(f']  {s.vertices[0].f: .6f} ')

print(__name__ + " module loaded", file=stderr)
