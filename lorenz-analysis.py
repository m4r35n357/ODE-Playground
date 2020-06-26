#!/usr/bin/env python3

#  Example: ./lorenz-analysis.py 3 -10 -10 10 -25 25 10.0 28.0 8.0 3.0
from sys import argv
from dual import Context
from ode_analysis import equilibrium, plot_lambda

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
    if len(argv) != 11:
        raise Exception(">>> ERROR! Please supply output precision, initial values, lambda range, and ODE parameters <<<")
    Context.places = argv[1]  # controls
    x, y, z = float(argv[2]), float(argv[3]), float(argv[4])  # initial values
    λa, λb = float(argv[5]), float(argv[6])
    σ, ρ, β = float(argv[7]), float(argv[8]), float(argv[9]) / float(argv[10])  # parameters
    print(equilibrium(f_a, f_b, f_c, x, y, z))
    plot_lambda(f_a, f_b, f_c, x, y, z, λ_min=λa, λ_max=λb, ce_min=-5000, ce_max=5000)

main()
