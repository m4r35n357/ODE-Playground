from sys import argv  # python3 lorenz.py 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3
order, δt, n_steps = int(argv[1]), float(argv[2]), int(argv[3])  # integrator controls
x, y, z = float(argv[4]), float(argv[5]), float(argv[6])  # coordinates
jx, jy, jz = [0.0] * (order + 1), [0.0] * (order + 1), [0.0] * (order + 1)  # jets
σ, ρ, β = float(argv[7]), float(argv[8]), float(argv[9]) / float(argv[10])  # parameters
for step in range(1, n_steps + 1):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * δt))
    jx[0], jy[0], jz[0] = x, y, z
    for k in range(order):  # build up the jets using the recurrence relations and the derivative rule
        jx[k + 1] = σ * (jy[k] - jx[k]) / (k + 1)
        jy[k + 1] = (ρ * jx[k] - sum(jx[j] * jz[k - j] for j in range(k + 1)) - jy[k]) / (k + 1)
        jz[k + 1] = (sum(jx[j] * jy[k - j] for j in range(k + 1)) - β * jz[k]) / (k + 1)
    x, y, z = jx[order], jy[order], jz[order]
    for i in range(order - 1, -1, -1):  # Horner's method
        x, y, z = x * δt + jx[i], y * δt + jy[i], z * δt + jz[i]
