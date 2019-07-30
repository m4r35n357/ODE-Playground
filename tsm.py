#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv
from math import sqrt, cos, pi
from ad import t_jet, t_horner, t_abs, t_prod, t_sin_cos, t_tan_sec2


def output(x, y, z, t):
    print(f'{x:.9e} {y:.9e} {z:.9e} {t:.5e}')

def main():
    model, order, δt, n_steps = argv[1], int(argv[2]), float(argv[3]), int(argv[4])  # integrator controls
    x0, y0, z0 = float(argv[5]), float(argv[6]), float(argv[7])  # initial values
    x, y, z = t_jet(order + 1), t_jet(order + 1), t_jet(order + 1)  # coordinate jets

    if model == "lorenz":
        #  Example: ./tsm.py lorenz 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        #  Example: ./tsm.py lorenz 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -25 50
        σ, ρ, β = float(argv[8]), float(argv[9]), float(argv[10]) / float(argv[11])
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = σ * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (ρ * x[k] - t_prod(x, z, k) - y[k]) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - β * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "lu":
        #  Example: ./tsm.py lu 10 .01 3000 -3 2 20 36 3 20 | ./plotPi3d.py
        #  Example: ./tsm.py lu 10 .01 3000 -3 2 20 36 3 20 | ./plotAnimated.py 1 -25 50
        a, b, c = float(argv[8]), float(argv[9]), float(argv[10])
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "chen":
        #  Example: ./tsm.py chen 10 .01 3000 -3 2 20 35 3 28 | ./plotPi3d.py
        #  Example: ./tsm.py chen 10 .01 3000 -3 2 20 35 3 28 | ./plotAnimated.py 1 -25 50
        a, b, c = float(argv[8]), float(argv[9]), float(argv[10])
        d = c - a
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (d * x[k] + c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rossler":
        #  Example: ./tsm.py rossler 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotPi3d.py
        #  Example: ./tsm.py rossler 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotAnimated.py 1 -20 30
        a, b_, c = float(argv[8]), t_jet(order, float(argv[9])), float(argv[10])
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = - (y[k] + z[k]) / (k + 1)
                y[k + 1] = (x[k] + a * y[k]) / (k + 1)
                z[k + 1] = (b_[k] + t_prod(x, z, k) - c * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "bouali":
        #  Example: ./tsm.py bouali 10 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotPi3d.py
        #  Example: ./tsm.py bouali 10 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotAnimated.py 1 -5 5
        a, b, c, d = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
        jet1, w4, w5 = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                w4[k] = jet1[k] - y[k]
                w5[k] = jet1[k] - t_prod(x, x, k)
                x[k + 1] = (a * t_prod(x, w4, k) - b * z[k]) / (k + 1)
                y[k + 1] = - c * t_prod(y, w5, k) / (k + 1)
                z[k + 1] = d * x[k] / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "thomas":
        #  Example: ./tsm.py thomas 10 0.1 30000 1 0 0 .19 | ./plotPi3d.py
        #  Example: ./tsm.py thomas 10 0.1 30000 1 0 0 .19 | ./plotAnimated.py 1 -5 5
        b = float(argv[8])
        sx, cx = t_jet(order), t_jet(order)
        sy, cy = t_jet(order), t_jet(order)
        sz, cz = t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                sy[k], cy[k] = t_sin_cos(sy, cy, y, k)
                sz[k], cz[k] = t_sin_cos(sz, cz, z, k)
                x[k + 1] = (sy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (sz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (sx[k] - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott-thomas":
        #  Example: ./tsm.py sprott-thomas 10 0.02 30000 1 0 0 4.75 .7 | ./plotPi3d.py
        #  Example: ./tsm.py sprott-thomas 10 0.02 30000 1 0 0 4.75 .7 | ./plotAnimated.py 1 -1 1
        a, b = float(argv[8]), float(argv[9])
        sx, cx, ax, sax, cax = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        sy, cy, ay, say, cay = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        sz, cz, az, saz, caz = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                sx[k], cx[k] = t_tan_sec2(sx, cx, x, k)
                sy[k], cy[k] = t_tan_sec2(sy, cy, y, k)
                sz[k], cz[k] = t_tan_sec2(sz, cz, z, k)
                ax[k], ay[k], az[k] = a * x[k], a * y[k], a * z[k]
                sax[k], cax[k] = t_sin_cos(sax, cax, ax, k)
                say[k], cay[k] = t_sin_cos(say, cay, ay, k)
                saz[k], caz[k] = t_sin_cos(saz, caz, az, k)
                x[k + 1] = (say[k] - b * sx[k]) / (k + 1)
                y[k + 1] = (saz[k] - b * sy[k]) / (k + 1)
                z[k + 1] = (sax[k] - b * sz[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rabinovich–fabrikant":
        #  Example: ./tsm.py rabinovich–fabrikant 10 .01 100001 .1 .1 .1 .2876 .1 | ./plotPi3d.py
        #  Example: ./tsm.py rabinovich–fabrikant 10 .01 100001 .1 .1 .1 .2876 .1 | ./plotAnimated.py 1 -3 3
        α, γ = float(argv[8]), float(argv[9])
        jet1, a, b, c = t_jet(order, 1), t_jet(order), t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x2_1 = t_prod(x, x, k) - jet1[k]
                a[k] = z[k] + x2_1
                b[k] = 3.0 * z[k] - x2_1
                c[k] = α + t_prod(x, y, k)
                x[k + 1] = (t_prod(y, a, k) + γ * x[k]) / (k + 1)
                y[k + 1] = (t_prod(x, b, k) + γ * y[k]) / (k + 1)
                z[k + 1] = - 2.0 * t_prod(z, c, k) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott":
        #  Example: ./tsm.py sprott 10 0.1 30001 1 0 0 | ./plotPi3d.py
        w1 = t_jet(order, 1)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = (y[k] + 2.0 * t_prod(x, y, k) + t_prod(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - 2.0 * t_prod(x, x, k) + t_prod(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - t_prod(x, x, k) - t_prod(y, y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott-jafari":
        #  Example: ./tsm.py sprott-jafari 10 0.01 30001 0 3.9 .7 8.888 4 | ./plotPi3d.py
        #  Example: ./tsm.py sprott-jafari 10 0.01 30001 0 3.9 .7 8.888 4 | ./plotAnimated.py 1 -20 20
        a, b_ = float(argv[8]), t_jet(order, float(argv[9]))
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - x[k] + t_prod(y, z, k) / (k + 1)
                z[k + 1] = (z[k] + a * t_prod(x, x, k) - t_prod(y, y, k) - b_[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "halvorsen":
        #  Example: ./tsm.py halvorsen 10 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        #  Example: ./tsm.py halvorsen 10 .01 100001 1 0 0 1.4 | ./plotAnimated.py 1 -15 10
        α = float(argv[8])
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = - (α * x[k] + 4.0 * y[k] + 4.0 * z[k] + t_prod(y, y, k)) / (k + 1)
                y[k + 1] = - (α * y[k] + 4.0 * z[k] + 4.0 * x[k] + t_prod(z, z, k)) / (k + 1)
                z[k + 1] = - (α * z[k] + 4.0 * x[k] + 4.0 * y[k] + t_prod(x, x, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "nose-hoover":
        #  Example: ./tsm.py nose-hoover 10 0.01 10001 1 0 0 6.0 | ./plotPi3d.py
        #  Example: ./tsm.py nose-hoover 10 0.01 10001 1 0 0 6.0 | ./plotAnimated.py 1 -10 10
        α = t_jet(order, float(argv[8]))
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (t_prod(y, z, k) - x[k]) / (k + 1)
                z[k + 1] = (α[k] - t_prod(y, y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rucklidge":
        #  Example: ./tsm.py rucklidge 10 0.01 10001 1 0 0 6.7 2 | ./plotPi3d.py
        #  Example: ./tsm.py rucklidge 10 0.01 10001 1 0 0 6.7 2 | ./plotAnimated.py 1 -15 20
        α, κ = float(argv[8]), float(argv[9])
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = (α * y[k] - κ * x[k] - t_prod(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_prod(y, y, k) - z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "wimol-banlue":
        #  Example: ./tsm.py wimol-banlue 10 0.1 10001 1 0 0 2.0 | ./plotPi3d.py
        #  Example: ./tsm.py wimol-banlue 10 0.1 10001 1 0 0 2.0 | ./plotAnimated.py 1 -5 5
        α = t_jet(order, float(argv[8]))
        tx, sx = t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                tx[k], sx[k] = t_tan_sec2(tx, sx, x, k, hyp=True)
                x[k + 1] = (y[k] - x[k]) / (k + 1)
                y[k + 1] = - t_prod(z, tx, k) / (k + 1)
                z[k + 1] = (- α[k] + t_prod(x, y, k) + t_abs(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "newton":
        κ, l, m = float(argv[8]), float(argv[9]), float(argv[10])
        output(x0, y0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (- x[k] - κ * m / l**2) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "oscillator":
        #  Example: ./tsm.py oscillator 10 .05 4001 0.0 0.0 0.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        κ, ζ, a, ω = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
        output(x0, y0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * y[k] - κ * x[k]) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "pendulum":
        #  Example: ./tsm.py pendulum 10 .05 4001 0.0 0.0 0.0 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        g, m, length = 9.80665, float(argv[8]), float(argv[9])  # physical parameters
        ζ, a, ω = float(argv[10]), float(argv[11]), 2.0 * pi * sqrt(length / g) * float(argv[12])  # damping/forcing
        sinθ, cosθ = t_jet(order), t_jet(order)  # jets
        output(x0, y0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):  # build up jets using recurrences and the derivative rule
                sinθ[k], cosθ[k] = t_sin_cos(sinθ, cosθ, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * length * y[k] - m * g * sinθ[k]) / (m * length) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)  # Horner's method
            output(x0, y0, 0.0, step * δt)
    elif model == "volterra":
        #  Example: ./tsm.py volterra 10 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 1 0 80
        a, b, c, d = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
        output(x0, y0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                xy = t_prod(x, y, k)
                x[k + 1] = (a * x[k] - c * xy) / (k + 1)
                y[k + 1] = (d * xy - b * y[k]) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "logistic":
        #  Example: ./tsm.py logistic 10 0.1 10001 .6 0 0 .1 | ./plotXY.py 1 3 0
        a = float(argv[8])
        w1, wa, wb = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, 0.0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0] = x0
            for k in range(order):
                wa[k] = a * x[k]
                wb[k] = w1[k] - x[k]
                x[k + 1] = t_prod(wa, wb, k) / (k + 1)
            x0 = t_horner(x, δt)
            output(x0, 0.0, 0.0, step * δt)
    elif model == "constant":
        #  Example: ./tsm.py constant 10 0.1 10001 10 0 0 -.05 | ./plotXY.py 1 3 0
        #  Example: ./tsm.py constant 10 0.1 10001 10 0 0 -.05 | ./plotAnimated.py 1 0 10
        a = float(argv[8])
        output(x0, 0.0, 0.0, 0.0)
        for step in range(1, n_steps + 1):
            x[0] = x0
            for k in range(order):
                x[k + 1] = a * x[k] / (k + 1)
            x0 = t_horner(x, δt)
            output(x0, 0.0, 0.0, step * δt)

main()
