#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv
from math import sqrt, cos, pi
from ad import Context, t_jet, t_horner, t_abs, t_prod, t_sin_cos, t_tan_sec2, t_sqr


def output(x, y, z, t):
    print(f'{x:+.{Context.places}e} {y:+.{Context.places}e} {z:+.{Context.places}e} {t:.5e}')

def main():
    model, Context.places, order, δt, n_steps = argv[1], int(argv[2]), int(argv[3]), float(argv[4]), int(argv[5])  # controls
    x0, y0, z0 = float(argv[6]), float(argv[7]), float(argv[8])  # initial values
    x, y, z = t_jet(order + 1), t_jet(order + 1), t_jet(order + 1)  # coordinate jets
    steps = range(1, n_steps + 1)
    index = range(order)

    if model == "lorenz":
        #  Example: ./tsm.py lorenz 9 8 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        #  Example: ./tsm.py lorenz 9 8 .01 10000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -25 50
        #  Profile: py-spy --function -- python tsm.py lorenz 9 8 .01 10000000 -15.8 -17.48 35.64 10 28 8 3
        σ, ρ, β = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = σ * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (ρ * x[k] - t_prod(x, z, k) - y[k]) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - β * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "lu":
        #  Example: ./tsm.py lu 9 8 .01 3000 -3 2 20 36 3 20 | ./plotPi3d.py
        #  Example: ./tsm.py lu 9 8 .01 3000 -3 2 20 36 3 20 | ./plotAnimated.py 1 -25 50
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        output(x0, y0, z0, 0.0)
        for step in steps:
            x, y, z = t_jet(order + 1), t_jet(order + 1), t_jet(order + 1)  # reset coordinate jets (for more visual debugging)
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "chen":
        #  Example: ./tsm.py chen 9 8 .01 3000 -3 2 20 35 3 28 | ./plotPi3d.py
        #  Example: ./tsm.py chen 9 8 .01 3000 -3 2 20 35 3 28 | ./plotAnimated.py 1 -25 50
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        d = c - a
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (d * x[k] + c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rossler":
        #  Example: ./tsm.py rossler 9 8 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotPi3d.py
        #  Example: ./tsm.py rossler 9 8 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotAnimated.py 1 -20 30
        a, b, c = float(argv[9]), t_jet(order, float(argv[10])), float(argv[11])
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = - (y[k] + z[k]) / (k + 1)
                y[k + 1] = (x[k] + a * y[k]) / (k + 1)
                z[k + 1] = (b[k] + t_prod(x, z, k) - c * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "bouali":
        #  Example: ./tsm.py bouali 9 8 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotPi3d.py
        #  Example: ./tsm.py bouali 9 8 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotAnimated.py 1 -5 5
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        jet1, w4, w5 = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                w4[k] = jet1[k] - y[k]
                w5[k] = jet1[k] - t_sqr(x, k)
                x[k + 1] = (a * t_prod(x, w4, k) - b * z[k]) / (k + 1)
                y[k + 1] = - c * t_prod(y, w5, k) / (k + 1)
                z[k + 1] = d * x[k] / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "thomas":
        #  Example: ./tsm.py thomas 9 8 0.1 30000 1 0 0 .19 | ./plotPi3d.py
        #  Example: ./tsm.py thomas 9 8 0.1 30000 1 0 0 .19 | ./plotAnimated.py 1 -5 5
        b = float(argv[9])
        sx, cx = t_jet(order), t_jet(order)
        sy, cy = t_jet(order), t_jet(order)
        sz, cz = t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                sy[k], cy[k] = t_sin_cos(sy, cy, y, k)
                sz[k], cz[k] = t_sin_cos(sz, cz, z, k)
                x[k + 1] = (sy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (sz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (sx[k] - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott-thomas-1" or model == "sprott-thomas-2":
        #  Example: ./tsm.py sprott-thomas-1 9 8 0.02 30000 1 0 0 4.75 1 | ./plotPi3d.py
        #  Example: ./tsm.py sprott-thomas-1 9 8 0.02 30000 1 0 0 4.75 1 | ./plotAnimated.py 1 -1 1
        #  Example: ./tsm.py sprott-thomas-2 9 8 0.02 30000 1 0 0 4.75 .7 | ./plotPi3d.py
        #  Example: ./tsm.py sprott-thomas-2 9 8 0.02 30000 1 0 0 4.75 .7 | ./plotAnimated.py 1 -1 1
        a, b = float(argv[9]), float(argv[10])
        sx, cx, ax, sax, cax = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        sy, cy, ay, say, cay = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        sz, cz, az, saz, caz = t_jet(order), t_jet(order), t_jet(order), t_jet(order), t_jet(order)
        fun = t_sin_cos if model == "sprott-thomas-1" else t_tan_sec2
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
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
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rabinovich–fabrikant":
        #  Example: ./tsm.py rabinovich–fabrikant 9 8 .01 100001 .1 .1 .1 .2876 .1 | ./plotPi3d.py
        #  Example: ./tsm.py rabinovich–fabrikant 9 8 .01 100001 .1 .1 .1 .2876 .1 | ./plotAnimated.py 1 -3 3
        α, γ = t_jet(order, float(argv[9])), float(argv[10])
        jet1, a, b, c = t_jet(order, 1.0), t_jet(order), t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x2_1 = t_sqr(x, k) - jet1[k]
                a[k] = z[k] + x2_1
                b[k] = 3.0 * z[k] - x2_1
                c[k] = α[k] + t_prod(x, y, k)
                x[k + 1] = (t_prod(y, a, k) + γ * x[k]) / (k + 1)
                y[k + 1] = (t_prod(x, b, k) + γ * y[k]) / (k + 1)
                z[k + 1] = - 2.0 * t_prod(z, c, k) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott":
        #  Example: ./tsm.py sprott 9 8 0.1 30001 1 0 0 | ./plotPi3d.py
        #  Example: ./tsm.py sprott 9 8 0.1 30001 1 0 0 | ./plotAnimated.py 1 -20 20
        w1 = t_jet(order, 1)
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x2 = t_sqr(x, k)
                x[k + 1] = (y[k] + 2.0 * t_prod(x, y, k) + t_prod(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - 2.0 * x2 + t_prod(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - x2 - t_sqr(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "sprott-jafari":
        #  Example: ./tsm.py sprott-jafari 9 8 0.01 30001 0 3.9 .7 8.888 4 | ./plotPi3d.py
        #  Example: ./tsm.py sprott-jafari 9 8 0.01 30001 0 3.9 .7 8.888 4 | ./plotAnimated.py 1 -20 20
        a, b = float(argv[9]), t_jet(order, float(argv[10]))
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - x[k] + t_prod(y, z, k) / (k + 1)
                z[k + 1] = (z[k] + a * t_sqr(x, k) - t_sqr(y, k) - b[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "halvorsen":
        #  Example: ./tsm.py halvorsen 9 8 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        #  Example: ./tsm.py halvorsen 9 8 .01 100001 1 0 0 1.4 | ./plotAnimated.py 1 -15 10
        α = float(argv[9])
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = - (α * x[k] + 4.0 * y[k] + 4.0 * z[k] + t_sqr(y, k)) / (k + 1)
                y[k + 1] = - (α * y[k] + 4.0 * z[k] + 4.0 * x[k] + t_sqr(z, k)) / (k + 1)
                z[k + 1] = - (α * z[k] + 4.0 * x[k] + 4.0 * y[k] + t_sqr(x, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "nose-hoover":
        #  Example: ./tsm.py nose-hoover 9 8 0.01 10001 1 0 0 6.0 | ./plotPi3d.py
        #  Example: ./tsm.py nose-hoover 9 8 0.01 10001 1 0 0 6.0 | ./plotAnimated.py 1 -10 10
        α = t_jet(order, float(argv[9]))
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (t_prod(y, z, k) - x[k]) / (k + 1)
                z[k + 1] = (α[k] - t_sqr(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "rucklidge":
        #  Example: ./tsm.py rucklidge 9 8 0.01 10001 1 0 0 6.7 2 | ./plotPi3d.py
        #  Example: ./tsm.py rucklidge 9 8 0.01 10001 1 0 0 6.7 2 | ./plotAnimated.py 1 -15 20
        α, κ = float(argv[9]), float(argv[10])
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                x[k + 1] = (α * y[k] - κ * x[k] - t_prod(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_sqr(y, k) - z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "wimol-banlue":
        #  Example: ./tsm.py wimol-banlue 9 8 0.1 10001 1 0 0 2.0 | ./plotPi3d.py
        #  Example: ./tsm.py wimol-banlue 9 8 0.1 10001 1 0 0 2.0 | ./plotAnimated.py 1 -5 5
        α = t_jet(order, float(argv[9]))
        tx, sx = t_jet(order), t_jet(order)
        output(x0, y0, z0, 0.0)
        for step in steps:
            x[0], y[0], z[0] = x0, y0, z0
            for k in index:
                tx[k], sx[k] = t_tan_sec2(tx, sx, x, k, hyp=True)
                x[k + 1] = (y[k] - x[k]) / (k + 1)
                y[k + 1] = - t_prod(z, tx, k) / (k + 1)
                z[k + 1] = (- α[k] + t_prod(x, y, k) + t_abs(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
            output(x0, y0, z0, step * δt)
    elif model == "oscillator":
        #  Example: ./tsm.py oscillator 9 8 .05 4001 0.0 0.0 0.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        κ, ζ, a, ω = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        output(x0, y0, 0.0, 0.0)
        for step in steps:
            x[0], y[0] = x0, y0
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * (2.0 * sqrt(κ)) * y[k] - κ * x[k]) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "pendulum":
        #  Example: ./tsm.py pendulum 9 8 .05 4001 0.0 0.0 0.0 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        g, m, length = 9.80665, float(argv[9]), float(argv[10])  # physical parameters
        ζ, a, ω = float(argv[11]), float(argv[12]), 2.0 * pi * sqrt(length / g) * float(argv[12])  # damping/forcing
        sinθ, cosθ = t_jet(order), t_jet(order)  # jets
        output(x0, y0, 0.0, 0.0)
        for step in steps:
            x[0], y[0] = x0, y0
            for k in index:  # build up jets using recurrences and the derivative rule
                sinθ[k], cosθ[k] = t_sin_cos(sinθ, cosθ, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * length * y[k] - m * g * sinθ[k]) / (m * length) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)  # Horner's method
            output(x0, y0, 0.0, step * δt)
    elif model == "volterra":
        #  Example: ./tsm.py volterra 9 8 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 1 0 80
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        output(x0, y0, 0.0, 0.0)
        for step in steps:
            x[0], y[0] = x0, y0
            for k in index:
                xy = t_prod(x, y, k)
                x[k + 1] = (a * x[k] - c * xy) / (k + 1)
                y[k + 1] = (d * xy - b * y[k]) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "logistic":
        #  Example: ./tsm.py logistic 9 8 0.1 10001 .6 0 0 .1 | ./plotXY.py 1 3 0
        #  Example: ./tsm.py logistic 9 8 0.1 10001 .6 0 0 .1 | ./plotAnimated.py 1 0 2
        a = float(argv[9])
        w1, wa, wb = t_jet(order, 1), t_jet(order), t_jet(order)
        output(x0, 0.0, 0.0, 0.0)
        for step in steps:
            x[0] = x0
            for k in index:
                wa[k] = a * x[k]
                wb[k] = w1[k] - x[k]
                x[k + 1] = t_prod(wa, wb, k) / (k + 1)
            x0 = t_horner(x, δt)
            output(x0, 0.0, 0.0, step * δt)
    elif model == "damped":
        #  Example: ./tsm.py damped 9 8 .1 1001 10 0 0 1 .5 | ./plotXY.py 1 3 0
        #  Example: ./tsm.py damped 9 8 .1 1001 10 0 0 1 .5 | ./plotAnimated.py 1 -10 10
        κ, ζ = float(argv[9]), float(argv[10])
        output(x0, y0, 0.0, 0.0)
        for step in steps:
            x[0], y[0] = x0, y0
            for k in index:
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - (κ * x[k] + ζ * (2.0 * sqrt(κ)) * y[k]) / (k + 1)
            x0, y0 = t_horner(x, δt), t_horner(y, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "constant":
        #  Example: ./tsm.py constant 9 8 0.1 10001 10 0 0 -.05 | ./plotXY.py 1 3 0
        #  Example: ./tsm.py constant 9 8 0.1 10001 10 0 0 -.05 | ./plotAnimated.py 1 0 10
        a = float(argv[9])
        output(x0, 0.0, 0.0, 0.0)
        for step in steps:
            x[0] = x0
            for k in index:
                x[k + 1] = a * x[k] / (k + 1)
            x0 = t_horner(x, δt)
            output(x0, 0.0, 0.0, step * δt)
    else:
        raise RuntimeError(f"Unknown model: {model}")


main()
