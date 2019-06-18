#!/usr/bin/env python3

#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from gmpy2 import get_context, log, sqrt, cos
from math import pi
get_context().precision = int(int(argv[2]) * log(10.0) / log(2.0))  # Set this BEFORE importing any AD stuff!
from ad import t_jet, t_horner, t_prod, t_sin_cos, t_tan_sec2, to_mpfr


def output(x, y, z, t):
    print(f'{x:.9e} {y:.9e} {z:.9e} {t:.5e}')


def main():
    d0 = to_mpfr(0)
    model, order, δt, n_steps = argv[1], int(argv[3]), to_mpfr(argv[4]), int(argv[5])  # integrator controls
    x0, y0, z0 = to_mpfr(argv[6]), to_mpfr(argv[7]), to_mpfr(argv[8])  # initial values
    x, y, z = t_jet(order + 1), t_jet(order + 1), t_jet(order + 1)  # coordinate jets

    if model == "lorenz":
        #  Example: ./tsm.py lorenz 16 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        #  Example: ./tsm.py lorenz 16 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -25 50
        σ, ρ, β = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]) / to_mpfr(argv[12])
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = σ * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (ρ * x[k] - t_prod(x, z, k) - y[k]) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - β * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "lu":
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "chen":
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        d = c - a
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (d * x[k] + c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rossler":
        #  Example: ./tsm.py rossler 16 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotPi3d.py
        #  Example: ./tsm.py rossler 16 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotAnimated.py 1 -20 30
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        b_ = t_jet(order, b)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = - (y[k] + z[k]) / (k + 1)
                y[k + 1] = (x[k] + a * y[k]) / (k + 1)
                z[k + 1] = (b_[k] + t_prod(x, z, k) - c * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "bouali":
        #  Example: ./tsm.py bouali 80 40 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotPi3d.py
        a, b, c, d = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]), to_mpfr(argv[12])
        jet1, w4, w5 = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                w4[k] = jet1[k] - y[k]
                w5[k] = jet1[k] - t_prod(x, x, k)
                x[k + 1] = (a * t_prod(x, w4, k) - b * z[k]) / (k + 1)
                y[k + 1] = - c * t_prod(y, w5, k) / (k + 1)
                z[k + 1] = d * x[k] / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "thomas":
        #  Example: ./tsm.py thomas 16 10 0.1 30000 1 0 0 .19 | ./plotPi3d.py
        #  Example: ./tsm.py thomas 16 10 0.1 30000 1 0 0 .19 | ./plotAnimated.py 1 -5 5
        b = to_mpfr(argv[9])
        sx, cx = t_jet(order), t_jet(order)
        sy, cy = t_jet(order), t_jet(order)
        sz, cz = t_jet(order), t_jet(order)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                sy[k], cy[k] = t_sin_cos(sy, cy, y, k)
                sz[k], cz[k] = t_sin_cos(sz, cz, z, k)
                x[k + 1] = (sy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (sz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (sx[k] - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "st":
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        sx, cx = t_jet(order), t_jet(order)
        sy, cy = t_jet(order), t_jet(order)
        sz, cz = t_jet(order), t_jet(order)
        ax, ay, az = t_jet(order), t_jet(order), t_jet(order)
        sax, cax = t_jet(order), t_jet(order)
        say, cay = t_jet(order), t_jet(order)
        saz, caz = t_jet(order), t_jet(order)
        output(x0, y0, z0, d0)
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
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rf":
        #  Example: ./tsm.py rf 16 10 .01 100001 .1 .1 .1 .2876 .1 | ./plotPi3d.py
        a, g = to_mpfr(argv[9]), to_mpfr(argv[10])
        jet1, w_a, w_b, w_c = t_jet(order, 1), t_jet(order), t_jet(order), t_jet(order)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                w_x2_1 = t_prod(x, x, k) - jet1[k]
                w_a[k] = z[k] + w_x2_1
                w_b[k] = 3 * z[k] - w_x2_1
                w_c[k] = a + t_prod(x, y, k)
                x[k + 1] = (t_prod(y, w_a, k) + g * x[k]) / (k + 1)
                y[k + 1] = (t_prod(x, w_b, k) + g * y[k]) / (k + 1)
                z[k + 1] = - 2 * t_prod(z, w_c, k) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott":
        #  Example: ./tsm.py sprott 16 10 0.1 30001 1 0 0 | ./plotPi3d.py
        w1 = t_jet(order, 1)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = (y[k] + 2 * t_prod(x, y, k) + t_prod(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - 2 * t_prod(x, x, k) + t_prod(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - t_prod(x, x, k) - t_prod(y, y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sj":
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        w_b = t_jet(order, b)
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - x[k] + t_prod(y, z, k) / (k + 1)
                z[k + 1] = (z[k] + a * t_prod(x, x, k) - t_prod(y, y, k) - w_b[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "halvorsen":
        #  Example: ./tsm.py halvorsen 16 10 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        a = to_mpfr(argv[9])
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = - (a * x[k] + 4 * y[k] + 4 * z[k] + t_prod(y, y, k)) / (k + 1)
                y[k + 1] = - (a * y[k] + 4 * z[k] + 4 * x[k] + t_prod(z, z, k)) / (k + 1)
                z[k + 1] = - (a * z[k] + 4 * x[k] + 4 * y[k] + t_prod(x, x, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "nh":
        a_ = t_jet(order, to_mpfr(argv[9]))
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (t_prod(y, z, k) - x[k]) / (k + 1)
                z[k + 1] = (a_[k] - t_prod(y, y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rucklidge":
        #  Example: ./tsm.py rucklidge 16 10 0.01 10001 1 0 0 6.7 2 | ./plotPi3d.py
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = (a * y[k] - b * x[k] - t_prod(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_prod(y, y, k) - z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "newton":
        κ, l, m = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        output(x0, y0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (- x[k] - κ * m / l**2) / (k + 1)
            x0, y0 = t_horner(x, order, δt), t_horner(y, order, δt)
            output(x0, y0, d0, step * δt)
    elif model == "oscillator":
            #  Example: ./tsm.py oscillator 16 10 .05 4001 0.0 0.0 0.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
            κ, ζ, a, ω = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]), to_mpfr(argv[12])
            output(x0, y0, d0, d0)
            for step in range(1, n_steps + 1):
                x[0], y[0] = x0, y0
                for k in range(order):
                    x[k + 1] = y[k] / (k + 1)
                    y[k + 1] = (a * cos(ω * step * δt) - ζ * y[k] - κ * x[k]) / (k + 1)
                x0, y0 = t_horner(x, order, δt), t_horner(y, order, δt)
                output(x0, y0, d0, step * δt)
    elif model == "pendulum":
        #  Example: ./tsm.py pendulum 16 10 .05 4001 0.0 0.0 0.0 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        g, m, length = 9.80665, to_mpfr(argv[9]), to_mpfr(argv[10])  # physical parameters
        ζ, a, ω = to_mpfr(argv[11]), to_mpfr(argv[12]), 2.0 * pi * sqrt(length / g) * to_mpfr(argv[13])  # damping/forcing
        sinθ, cosθ = t_jet(order), t_jet(order)  # jets
        output(x0, y0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):  # build up jets using recurrences and the derivative rule
                sinθ[k], cosθ[k] = t_sin_cos(sinθ, cosθ, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * length * y[k] - m * g * sinθ[k]) / (m * length) / (k + 1)
            x0, y0 = t_horner(x, order, δt), t_horner(y, order, δt)  # Horner's method
            output(x0, y0, d0, step * δt)
    elif model == "volterra":
        #  Example: ./tsm.py volterra 16 10 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 1 0 80
        a, b, c, d = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]), to_mpfr(argv[12])
        output(x0, y0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                xy = t_prod(x, y, k)
                x[k + 1] = (a * x[k] - c * xy) / (k + 1)
                y[k + 1] = (d * xy - b * y[k]) / (k + 1)
            x0, y0 = t_horner(x, order, δt), t_horner(y, order, δt)
            output(x0, y0, d0, step * δt)
    elif model == "logistic":
        #  Example: ./tsm.py logistic 16 10 0.1 10001 .6 0 0 .1 | ./plotXY.py 1 3 0
        a = to_mpfr(argv[9])
        w1, wa, wb = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, d0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0] = x0
            for k in range(order):
                wa[k] = a * x[k]
                wb[k] = w1[k] - x[k]
                x[k + 1] = t_prod(wa, wb, k) / (k + 1)
            x0 = t_horner(x, order, δt)
            output(x0, d0, d0, step * δt)
    elif model == "constant":
        #  Example: ./tsm.py constant 16 10 0.1 10001 10 0 0 -.05 | ./plotXY.py 1 3 0
        a = to_mpfr(argv[9])
        output(x0, d0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0] = x0
            for k in range(order):
                x[k + 1] = a * x[k] / (k + 1)
            x0 = t_horner(x, order, δt)
            output(x0, d0, d0, step * δt)


main()
