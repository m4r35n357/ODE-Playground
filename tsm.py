#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv
from math import sqrt, cos, pi, floor, sin
from ad import Context, tsm_jet, horner, t_abs, t_mul, t_sin_cos, t_tan_sec2, t_sqr, t_exp, t_div


def output(x, y, z, t):
    print(f'{x:+.{Context.places}e} {y:+.{Context.places}e} {z:+.{Context.places}e} {t:.5e}')

def main():
    model, Context.places, n, δt, n_steps = argv[1], int(argv[2]), int(argv[3]), float(argv[4]), int(argv[5])  # controls
    x, y, z = tsm_jet(n + 1), tsm_jet(n + 1), tsm_jet(n + 1)  # coordinate jets
    x[0], y[0], z[0] = float(argv[6]), float(argv[7]), float(argv[8])  # initial values
    steps = range(1, n_steps + 1)
    index = range(n)

    if model == "lorenz":
        #  Example: ./tsm.py lorenz 9 8 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plot3d.py
        #  Example: ./tsm.py lorenz 9 8 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py -25 50
        #  Profile: py-spy --function -- python tsm.py lorenz 9 8 .01 10000000 -15.8 -17.48 35.64 10 28 8 3
        σ, ρ, β = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = σ * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (ρ * x[k] - y[k] - t_mul(x, z, k)) / (k + 1)
                z[k + 1] = (t_mul(x, y, k) - β * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "lu" or model == "chen":
        #  Example: ./tsm.py lu 9 8 .01 3000 -3 2 20 36 3 20 | ./plot3d.py
        #  Example: ./tsm.py chen 9 8 .01 3000 -3 2 20 35 3 28 | ./plotAnimated.py -25 50
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        d = c - a if model == "chen" else 0.0
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            x, y, z = tsm_jet(n + 1), tsm_jet(n + 1), tsm_jet(n + 1)  # reset coordinate jets (for more visual debugging)
            for k in index:
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (d * x[k] + c * y[k] - t_mul(x, z, k)) / (k + 1)
                z[k + 1] = (t_mul(x, y, k) - b * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "rossler":
        #  Example: ./tsm.py rossler 9 8 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plot3d.py
        #  Example: ./tsm.py rossler 9 8 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotAnimated.py -20 30
        a, b, c = float(argv[9]), tsm_jet(n, float(argv[10])), float(argv[11])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = - (y[k] + z[k]) / (k + 1)
                y[k + 1] = (x[k] + a * y[k]) / (k + 1)
                z[k + 1] = (b[k] + t_mul(x, z, k) - c * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "genesio-tesi":
        #  Example: ./tsm.py genesio-tesi 9 8 0.01 150000 .1 .1 .1 .44 1.1 1 | ./plot3d.py
        #  Example: ./tsm.py genesio-tesi 9 8 0.01 150000 0.0 -6.78 0.02 .44 1.1 1 | ./plotAnimated.py -20 30
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = z[k] / (k + 1)
                z[k + 1] = (t_sqr(x, k) - c * x[k] - b * y[k] - a * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "bouali":
        #  Example: ./tsm.py bouali 9 8 0.02 50001 1 1 0 3 2.2 1 .01 | ./plot3d.py
        #  Example: ./tsm.py bouali 9 8 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotAnimated.py -5 5
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        jet1, w4, w5 = tsm_jet(n, 1), tsm_jet(n), tsm_jet(n)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                w4[k] = jet1[k] - y[k]
                w5[k] = jet1[k] - t_sqr(x, k)
                x[k + 1] = (a * t_mul(x, w4, k) - b * z[k]) / (k + 1)
                y[k + 1] = - c * t_mul(y, w5, k) / (k + 1)
                z[k + 1] = d * x[k] / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "thomas":
        #  Example: ./tsm.py thomas 9 8 0.1 30000 1 0 0 .19 | ./plot3d.py
        #  Example: ./tsm.py thomas 9 8 0.1 30000 1 0 0 .19 | ./plotAnimated.py -5 5
        b = float(argv[9])
        sx, cx = tsm_jet(n), tsm_jet(n)
        sy, cy = tsm_jet(n), tsm_jet(n)
        sz, cz = tsm_jet(n), tsm_jet(n)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                sy[k], cy[k] = t_sin_cos(sy, cy, y, k)
                sz[k], cz[k] = t_sin_cos(sz, cz, z, k)
                x[k + 1] = (sy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (sz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (sx[k] - b * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "sprott-thomas-1" or model == "sprott-thomas-2":
        #  Example: ./tsm.py sprott-thomas-1 9 8 0.02 30000 1 0 0 4.75 1 | ./plot3d.py
        #  Example: ./tsm.py sprott-thomas-1 9 8 0.02 30000 1 0 0 4.75 1 | ./plotAnimated.py -1 1
        #  Example: ./tsm.py sprott-thomas-2 9 8 0.02 30000 1 0 0 4.75 .7 | ./plot3d.py
        #  Example: ./tsm.py sprott-thomas-2 9 8 0.02 30000 1 0 0 4.75 .7 | ./plotAnimated.py -1 1
        a, b = float(argv[9]), float(argv[10])
        sx, cx, ax, sax, cax = tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n)
        sy, cy, ay, say, cay = tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n)
        sz, cz, az, saz, caz = tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n)
        fun = t_sin_cos if model == "sprott-thomas-1" else t_tan_sec2
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                sx[k], cx[k] = fun(sx, cx, x, k)
                sy[k], cy[k] = fun(sy, cy, y, k)
                sz[k], cz[k] = fun(sz, cz, z, k)
                ax[k], ay[k], az[k] = a * x[k], a * y[k], a * z[k]
                sax[k], cax[k] = t_sin_cos(sax, cax, ax, k)
                say[k], cay[k] = t_sin_cos(say, cay, ay, k)
                saz[k], caz[k] = t_sin_cos(saz, caz, az, k)
                x[k + 1] = (say[k] - b * sx[k]) / (k + 1)
                y[k + 1] = (saz[k] - b * sy[k]) / (k + 1)
                z[k + 1] = (sax[k] - b * sz[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "rabinovich–fabrikant":
        #  Example: ./tsm.py rabinovich–fabrikant 9 8 .01 100001 .1 .1 .1 .2876 .1 | ./plot3d.py
        #  Example: ./tsm.py rabinovich–fabrikant 9 8 .01 100001 .1 .1 .1 .2876 .1 | ./plotAnimated.py -3 3
        α, γ = tsm_jet(n, float(argv[9])), float(argv[10])
        jet1, a, b, c = tsm_jet(n, 1.0), tsm_jet(n), tsm_jet(n), tsm_jet(n)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x2_1 = t_sqr(x, k) - jet1[k]
                a[k] = z[k] + x2_1
                b[k] = 3.0 * z[k] - x2_1
                c[k] = α[k] + t_mul(x, y, k)
                x[k + 1] = (t_mul(y, a, k) + γ * x[k]) / (k + 1)
                y[k + 1] = (t_mul(x, b, k) + γ * y[k]) / (k + 1)
                z[k + 1] = - 2.0 * t_mul(z, c, k) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "sprott":
        #  Example: ./tsm.py sprott 9 8 0.1 30001 1 0 0 | ./plot3d.py
        #  Example: ./tsm.py sprott 9 8 0.1 30001 1 0 0 | ./plotAnimated.py -20 20
        w1 = tsm_jet(n, 1)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x2 = t_sqr(x, k)
                x[k + 1] = (y[k] + 2.0 * t_mul(x, y, k) + t_mul(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - 2.0 * x2 + t_mul(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - x2 - t_sqr(y, k)) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "sprott-jafari":
        #  Example: ./tsm.py sprott-jafari 9 8 0.01 30001 0 3.9 .7 8.888 4 | ./plot3d.py
        #  Example: ./tsm.py sprott-jafari 9 8 0.01 30001 0 3.9 .7 8.888 4 | ./plotAnimated.py -20 20
        a, b = float(argv[9]), tsm_jet(n, float(argv[10]))
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - x[k] + t_mul(y, z, k) / (k + 1)
                z[k + 1] = (z[k] + a * t_sqr(x, k) - t_sqr(y, k) - b[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "sprott-minimal":
        #  Example: ./tsm.py sprott-minimal 9 8 0.01 10001 .02 0 0 2.017 | ./plot3d.py
        #  Example: ./tsm.py sprott-minimal 9 8 0.01 10001 .02 0 0 2.017 | ./plotAnimated.py -10 10
        α = float(argv[9])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = z[k] / (k + 1)
                z[k + 1] = (- α * z[k] + t_sqr(y, k) - x[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "halvorsen":
        #  Example: ./tsm.py halvorsen 9 8 .01 100001 1 0 0 1.4 | ./plot3d.py
        #  Example: ./tsm.py halvorsen 9 8 .01 100001 1 0 0 1.4 | ./plotAnimated.py -15 10
        α = float(argv[9])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = - (α * x[k] + 4.0 * y[k] + 4.0 * z[k] + t_sqr(y, k)) / (k + 1)
                y[k + 1] = - (α * y[k] + 4.0 * z[k] + 4.0 * x[k] + t_sqr(z, k)) / (k + 1)
                z[k + 1] = - (α * z[k] + 4.0 * x[k] + 4.0 * y[k] + t_sqr(x, k)) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "nose-hoover":
        #  Example: ./tsm.py nose-hoover 9 8 0.01 10001 1 0 0 6.0 | ./plot3d.py
        #  Example: ./tsm.py nose-hoover 9 8 0.01 10001 1 0 0 6.0 | ./plotAnimated.py -10 10
        α = tsm_jet(n, float(argv[9]))
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (t_mul(y, z, k) - x[k]) / (k + 1)
                z[k + 1] = (α[k] - t_sqr(y, k)) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "rucklidge":
        #  Example: ./tsm.py rucklidge 9 8 0.01 10001 1 0 0 6.7 2 | ./plot3d.py
        #  Example: ./tsm.py rucklidge 9 8 0.01 10001 1 0 0 6.7 2 | ./plotAnimated.py -15 20
        α, κ = float(argv[9]), float(argv[10])
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = (α * y[k] - κ * x[k] - t_mul(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_sqr(y, k) - z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "wimol-banlue":
        #  Example: ./tsm.py wimol-banlue 9 8 0.1 10001 1 0 0 2.0 | ./plot3d.py
        #  Example: ./tsm.py wimol-banlue 9 8 0.1 10001 1 0 0 2.0 | ./plotAnimated.py -5 5
        α = tsm_jet(n, float(argv[9]))
        tx, sx = tsm_jet(n), tsm_jet(n)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                tx[k], sx[k] = t_tan_sec2(tx, sx, x, k, hyp=True)
                x[k + 1] = (y[k] - x[k]) / (k + 1)
                y[k + 1] = - t_mul(z, tx, k) / (k + 1)
                z[k + 1] = (- α[k] + t_mul(x, y, k) + t_abs(y, k)) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "yu-wang":
        #  Example: ./tsm.py yu-wang 9 8 .001 50000 1 0 0 10 40 2 2.5 | ./plot3d.py
        #  Example: ./tsm.py yu-wang 9 8 .001 50000 1 0 0 10 40 2 2.5 | ./plotAnimated.py -25 50
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        xy, e_xy = tsm_jet(n), tsm_jet(n)
        output(x[0], y[0], z[0], 0.0)
        for step in steps:
            for k in index:
                xy[k] = t_mul(x, y, k)
                e_xy[k] = t_exp(e_xy, xy, k)
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (b * x[k] - c * t_mul(x, z, k)) / (k + 1)
                z[k + 1] = (e_xy[k] - d * z[k]) / (k + 1)
            x[0], y[0], z[0] = horner(x, δt), horner(y, δt), horner(z, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "oscillator":
        #  Example: ./tsm.py oscillator 9 8 .05 4001 0.0 0.0 0.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py -50 50
        κ, ζ, a, ω = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        output(x[0], y[0], 0.0, 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * (2.0 * sqrt(κ)) * y[k] - κ * x[k]) / (k + 1)
            x[0], y[0] = horner(x, δt), horner(y, δt)
            output(x[0], y[0], z[0], step * δt)
    elif model == "pendulum":
        #  Example: ./tsm.py pendulum 9 8 .05 4001 0.0 0.0 0.0 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py -50 50
        g, m, length = 9.80665, float(argv[9]), float(argv[10])  # physical parameters
        ζ, a, ω = float(argv[11]), float(argv[12]), 2.0 * pi * sqrt(length / g) * float(argv[12])  # damping/forcing
        sinθ, cosθ = tsm_jet(n), tsm_jet(n)  # jets
        output(x[0], y[0], 0.0, 0.0)
        for step in steps:
            for k in index:  # build up jets using recurrences and the derivative rule
                sinθ[k], cosθ[k] = t_sin_cos(sinθ, cosθ, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * length * y[k] - m * g * sinθ[k]) / (m * length) / (k + 1)
            x[0], y[0] = horner(x, δt), horner(y, δt)  # Horner's method
            output(x[0], y[0], z[0], step * δt)
    elif model == "double-pendulum":
        # Example: ./tsm.py double-pendulum 9 10 0.01 10000 1 1 1 1 3 -1 3 -1 | ./plotPi2d.py
        def polar_to_rectangular(l_1, l_2, t, th_1, th_2, w_1, w_2):
            x_1, y_1 = l_1 * sin(th_1), - l_1 * cos(th_1)  # convert angles to X-Y coordinates
            x_2, y_2 = x_1 + l_2 * sin(th_2), y_1 - l_2 * cos(th_2)
            print(f"{x_1:.9e} {y_1:.9e} {x_2:.9e} {y_2:.9e} {t:.5e} 0.0 {th_1:.9e} {th_2:.9e} {w_1:.9e} {w_2:.9e}")
        g, m1, m2, l1, l2, θ1_0, ω1_0, θ2_0, ω2_0 = 9.80665, float(argv[6]), float(argv[7]), float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12]), float(argv[13])  # physical parameters
        θ1, θ2, sinθ1, cosθ1 = tsm_jet(n + 1), tsm_jet(n + 1), tsm_jet(n), tsm_jet(n)  # jets
        _θ1_θ2, sinθ1_θ2, cosθ1_θ2 = tsm_jet(n), tsm_jet(n), tsm_jet(n)
        _θ1_2θ2, sinθ1_2θ2, cosθ1_2θ2 = tsm_jet(n), tsm_jet(n), tsm_jet(n)
        _2θ1_2θ2, sin2θ1_2θ2, cos2θ1_2θ2 = tsm_jet(n), tsm_jet(n), tsm_jet(n)
        ω1, ω2 = tsm_jet(n + 1), tsm_jet(n + 1)
        ω1_2, ω2_2 = tsm_jet(n), tsm_jet(n)
        n1, n2, n1_, n2_, q1, q2, d = tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n), tsm_jet(n)
        polar_to_rectangular(l1, l2, 0.0, θ1_0, θ2_0, ω1_0, ω2_0)
        for step in steps:
            θ1[0], θ2[0], ω1[0], ω2[0] = θ1_0, θ2_0, ω1_0, ω2_0
            for k in index:  # build up jets using recurrences and the derivative rule
                sinθ1[k], cosθ1[k] = t_sin_cos(sinθ1, cosθ1, θ1, k)
                _θ1_θ2[k] = θ1[k] - θ2[k]
                sinθ1_θ2[k], cosθ1_θ2[k] = t_sin_cos(sinθ1_θ2, cosθ1_θ2, _θ1_θ2, k)
                _θ1_2θ2[k] = θ1[k] - 2.0 * θ2[k]
                sinθ1_2θ2[k], cosθ1_2θ2[k] = t_sin_cos(sinθ1_2θ2, cosθ1_2θ2, _θ1_2θ2, k)
                _2θ1_2θ2[k] = 2.0 * (θ1[k] - θ2[k])
                sin2θ1_2θ2[k], cos2θ1_2θ2[k] = t_sin_cos(sin2θ1_2θ2, cos2θ1_2θ2, _2θ1_2θ2, k)
                ω1_2[k] = t_sqr(ω1, k)
                ω2_2[k] = t_sqr(ω2, k)
                n1_[k] = l2 * ω2_2[k] + l1 * t_mul(ω1_2, cosθ1_θ2, k)
                n1[k] = - g * (2.0 * m1 + m2) * sinθ1[k] - g * m2 * sinθ1_2θ2[k] - 2.0 * m2 * t_mul(sinθ1_θ2, n1_, k)
                n2_[k] = g * (m1 + m2) * cosθ1[k] + l1 * (m1 + m2) * ω1_2[k] + l2 * m2 * t_mul(ω2_2, cosθ1_θ2, k)
                n2[k] = 2.0 * t_mul(sinθ1_θ2, n2_, k)
                d[k] = 2.0 * m1 + m2 - m2 * cos2θ1_2θ2[k]
                θ1[k + 1] = ω1[k] / (k + 1)
                θ2[k + 1] = ω2[k] / (k + 1)
                ω1[k + 1] = t_div(q1, n1, d, k) / l1 / (k + 1)
                ω2[k + 1] = t_div(q2, n2, d, k) / l2 / (k + 1)
            θ1_0, θ2_0, ω1_0, ω2_0 = horner(θ1, δt), horner(θ2, δt), horner(ω1, δt), horner(ω2, δt)  # Horner's method
            polar_to_rectangular(l1, l2, step * δt, θ1_0, θ2_0, ω1_0, ω2_0)
    elif model == "volterra":
        #  Example: ./tsm.py volterra 9 8 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 0 80
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        output(x[0], y[0], 0.0, 0.0)
        for step in steps:
            for k in index:
                xy = t_mul(x, y, k)
                x[k + 1] = (a * x[k] - c * xy) / (k + 1)
                y[k + 1] = (d * xy - b * y[k]) / (k + 1)
            x[0], y[0] = horner(x, δt), horner(y, δt)
            if floor(x[0]) == 0 or floor(y[0]) == 0:
                break
            output(floor(x[0]), floor(y[0]), 0.0, step * δt)
    elif model == "logistic":
        #  Example: ./tsm.py logistic 9 8 0.1 10001 .6 0 0 .1 | ./plotXY.py 3 0
        #  Example: ./tsm.py logistic 9 8 0.1 10001 .6 0 0 .1 | ./plotAnimated.py 0 2
        a = float(argv[9])
        w1, wa, wb = tsm_jet(n, 1), tsm_jet(n), tsm_jet(n)
        output(x[0], 0.0, 0.0, 0.0)
        for step in steps:
            for k in index:
                wa[k] = a * x[k]
                wb[k] = w1[k] - x[k]
                x[k + 1] = t_mul(wa, wb, k) / (k + 1)
            x[0] = horner(x, δt)
            output(x[0], 0.0, 0.0, step * δt)
    elif model == "damped":
        #  Example: ./tsm.py damped 9 8 .1 1001 10 0 0 1 .5 | ./plotXY.py 3 0
        #  Example: ./tsm.py damped 9 8 .1 1001 10 0 0 1 .5 | ./plotAnimated.py -10 10
        κ, ζ = float(argv[9]), float(argv[10])
        output(x[0], y[0], 0.0, 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - (κ * x[k] + 2.0 * sqrt(κ) * ζ * y[k]) / (k + 1)
            x[0], y[0] = horner(x, δt), horner(y, δt)
            output(x[0], y[0], 0.0, step * δt)
    elif model == "constant":
        #  Example: ./tsm.py constant 9 8 0.1 10001 10 0 0 -.05 | ./plotXY.py 3 0
        #  Example: ./tsm.py constant 9 8 0.1 10001 10 0 0 -.05 | ./plotAnimated.py 0 10
        a = float(argv[9])
        output(x[0], 0.0, 0.0, 0.0)
        for step in steps:
            for k in index:
                x[k + 1] = a * x[k] / (k + 1)
            x[0] = horner(x, δt)
            output(x[0], 0.0, 0.0, step * δt)
    else:
        raise RuntimeError(f"Unknown model: {model}")


main()
