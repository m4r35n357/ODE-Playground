#!/usr/bin/env python3

from dual import  *

sigma, rho, beta = 10.0, 28.0, 8.0 / 3

def f_a(x, y, z):
    return sigma * (y - x)

def f_b(x, y, z):
    return x * (rho - z) - y

def f_c(x, y, z):
    return x * y - beta * z

print(equilibrium(f_a, f_b, f_c, 8.5, 8.5, 27.0)) 
