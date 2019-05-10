#!/usr/bin/env python3
#  Example: ./forced.py 10 .05 4001 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
from sys import argv
from math import pi, sqrt, sin, cos

def dd(v, u, i):
    return sum(j * u[j] * v[i - j] for j in range(1, i)) / i

order, δt, n_steps = int(argv[1]), float(argv[2]), int(argv[3])  # integrator controls
g, m, l = 9.80665, float(argv[4]), float(argv[5])  # physical parameters
ζ, a, ω = float(argv[6]), float(argv[7]), 2.0 * pi * sqrt(l / g) * float(argv[8])  # damping/forcing

θ0, θdot0 = 0.0, 0.0  # initial values
θ, θdot, sinθ, cosθ = [0.0] * (order + 1), [0.0] * (order + 1), [0.0] * order, [0.0] * order  # jets
for step in range(1, n_steps):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(θ0, θdot0, 0.0, step * δt))

    # build up jets using recurrences and the derivative rule
    θ[0], θdot[0] = θ0, θdot0
    for k in range(order):
        sinθ[k], cosθ[k] = (sin(θ[0]), cos(θ[0])) if k == 0 else (cosθ[0] * θ[k] + dd(cosθ, θ, k),
                                                                - sinθ[0] * θ[k] - dd(sinθ, θ, k))
        θ[k + 1] = θdot[k] / (k + 1)
        θdot[k + 1] = (a * cos(ω * step * δt) - ζ * l * θdot[k] - m * g * sinθ[k]) / (m * l) / (k + 1)

    # evaluate series using Horner's method
    θ0, θdot0 = θ[order], θdot[order]
    for k in range(order - 1, -1, -1):
        θ0, θdot0 = θ0 * δt + θ[k], θdot0 * δt + θdot[k]
