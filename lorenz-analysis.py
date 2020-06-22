#!/usr/bin/env python3

from sys import argv
from dual import Context, equilibrium

σ, ρ, β = 0.0, 0.0, 0.0

def f_a(x, y, z):
    return σ * (y - x)

def f_b(x, y, z):
    return x * (ρ - z) - y

def f_c(x, y, z):
    return x * y - β * z

def main():
    global σ, ρ, β
    print(f'Lorenz System Analysis: {argv}')
    if len(argv) != 9:
        raise Exception(">>> ERROR! Please supply output precision, three initial values, and ODE parameters <<<")
    Context.places = argv[1]  # controls
    x, y, z = float(argv[2]), float(argv[3]), float(argv[4])  # initial values
    σ, ρ, β = float(argv[5]), float(argv[6]), float(argv[7]) / float(argv[8])  # parameters
    print(equilibrium(f_a, f_b, f_c, x, y, z))

main()
