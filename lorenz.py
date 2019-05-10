#!/usr/bin/env python3

from sys import argv  # ./lorenz.py 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3

order = int(argv[1])  # integrator controls
δt = float(argv[2])
n_steps = int(argv[3])

x = float(argv[4])  # coordinates
y = float(argv[5])
z = float(argv[6])

σ = float(argv[7])  # parameters
ρ = float(argv[8])
β = float(argv[9]) / float(argv[10])

jx = [0.0] * (order + 1)  # jets
jy = [0.0] * (order + 1)
jz = [0.0] * (order + 1)

for step in range(1, n_steps + 1):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * δt))

    jx[0] = x
    jy[0] = y
    jz[0] = z
    for k in range(order):  # build up jets using recurrences and the derivative rule
        jx[k + 1] = σ * (jy[k] - jx[k]) / (k + 1)
        jy[k + 1] = (ρ * jx[k] - sum(jx[j] * jz[k - j] for j in range(k + 1)) - jy[k]) / (k + 1)
        jz[k + 1] = (sum(jx[j] * jy[k - j] for j in range(k + 1)) - β * jz[k]) / (k + 1)

    x = jx[order]
    y = jy[order]
    z = jz[order]
    for i in range(order - 1, -1, -1):  # Horner's method
        x = x * δt + jx[i]
        y = y * δt + jy[i]
        z = z * δt + jz[i]
