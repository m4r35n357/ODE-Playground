#!/usr/bin/env python3
#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from math import sqrt, exp, sinh, cosh, sin, cos, tanh, tan, log, pi


def t_jet(n, value=0.0):
    j = [0.0] * n
    j[0] = value
    return j


def t_horner(j, n, δt):
    result = j[n]
    for i in range(n - 1, -1, -1):
        result = result * δt + j[i]
    return result


def t_prod(u, v, k):
    return sum(u[j] * v[k - j] for j in range(k + 1))


def t_quot(q, u, v, k):
    return (u[k] - sum(v[j] * q[k - j] for j in range(1, k + 1))) / v[0]


def _dd(v, u, k):
    return sum(j * u[j] * v[k - j] for j in range(1, k)) / k


def t_sqrt(r, u, k):
    return sqrt(u[0]) if k == 0 else (u[k] / 2.0 - _dd(r, r, k)) / r[0]


def t_exp(e, u, k):
    return exp(u[0]) if k == 0 else e[0] * u[k] + _dd(e, u, k)


def t_sin_cos(s, c, u, k, hyperbolic=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyperbolic else (sin(u[0]), cos(u[0]))
    else:
        sn = c[0] * u[k] + _dd(c, u, k)
        cn = s[0] * u[k] + _dd(s, u, k)
        return sn, cn if hyperbolic else - cn


def t_tan_sec2(t, s2, u, k, hyperbolic=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyperbolic else (tan(u[0]), tan(u[0])**2 + 1.0)
    else:
        tn = s2[0] * u[k] + _dd(s2, u, k)
        s2 = 2.0 * (t[0] * tn + _dd(t, t, k))
        return tn, - s2 if hyperbolic else s2


def t_pwr(p, u, a, k):
    return u[0] ** a if k == 0 else (a * (p[0] * u[k] + _dd(p, u, k)) - _dd(u, p, k)) / u[0]


def t_ln(l, u, k):
    return log(u[0]) if k == 0 else (u[k] - _dd(u, l, k)) / u[0]


def output(x, y, z, t):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, t))


def main():
    model = argv[1]
    order = int(argv[3])
    δt = float(argv[4])
    n_steps = int(argv[5])
    d0 = 0.0

    x0, y0, z0 = float(argv[6]), float(argv[7]), float(argv[8])
    x, y, z = t_jet(order + 1), t_jet(order + 1), t_jet(order + 1)

    if model == "lorenz":
        #  Example: ./tsm-float.py lorenz 16 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        #  Example: ./tsm-float.py lorenz 16 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotAnimated.py 1 -25 50
        σ, ρ, β = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
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
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        #  Example: ./tsm-float.py rossler 16 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotPi3d.py
        #  Example: ./tsm-float.py rossler 16 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotAnimated.py 1 -20 30
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        #  Example: ./tsm-float.py bouali 80 40 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotPi3d.py
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
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
        #  Example: ./tsm-float.py thomas 16 10 0.1 30000 1 0 0 .19 | ./plotPi3d.py
        #  Example: ./tsm-float.py thomas 16 10 0.1 30000 1 0 0 .19 | ./plotAnimated.py 1 -5 5
        b = float(argv[9])
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
        a, b = float(argv[9]), float(argv[10])
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
        #  Example: ./tsm-float.py rf 16 10 .01 100001 .1 .1 .1 .2876 .1 | ./plotPi3d.py
        a, g = float(argv[9]), float(argv[10])
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
        #  Example: ./tsm-float.py sprott 16 10 0.1 30001 1 0 0 | ./plotPi3d.py
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
        a, b = float(argv[9]), float(argv[10])
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
        #  Example: ./tsm-float.py halvorsen 16 10 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        a = float(argv[9])
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
        a_ = t_jet(order, int(argv[9]))
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
        #  Example: ./tsm-float.py rucklidge 16 10 0.01 10001 1 0 0 6.7 2 | ./plotPi3d.py
        a, b = float(argv[9]), float(argv[10])
        output(x0, y0, z0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(order):
                x[k + 1] = (a * y[k] - b * x[k] - t_prod(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_prod(y, y, k) - z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, order, δt), t_horner(y, order, δt), t_horner(z, order, δt)
            output(x0, y0, z0, step * δt)
    elif model == "damped":
        #  Example: ./tsm-float.py damped 16 10 .05 4001 0.0 0.0 0.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        κ, ζ, a, ω = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        output(x0, y0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0], y[0] = x0, y0
            for k in range(order):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (a * cos(ω * step * δt) - ζ * y[k] - κ * x[k]) / (k + 1)
            x0, y0 = t_horner(x, order, δt), t_horner(y, order, δt)
            output(x0, y0, 0.0, step * δt)
    elif model == "forced":
        #  Example: ./tsm-float.py forced 16 10 .05 4001 0.0 0.0 0.0 1.0 1.0 0.1 4.9 1.1 | ./plotAnimated.py 1 -50 50
        g, m, length = 9.80665, float(argv[9]), float(argv[10])  # physical parameters
        ζ, a, ω = float(argv[11]), float(argv[12]), 2.0 * pi * sqrt(length / g) * float(argv[13])  # damping/forcing
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
        #  Example: ./tsm-float.py volterra 16 10 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 1 0 80
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
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
        #  Example: ./tsm-float.py logistic 16 10 0.1 10001 .6 0 0 .1 | ./plotXY.py 1 3 0
        a = float(argv[9])
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
        #  Example: ./tsm-float.py constant 16 10 0.1 10001 10 0 0 -.05 | ./plotXY.py 1 3 0
        a = float(argv[9])
        output(x0, d0, d0, d0)
        for step in range(1, n_steps + 1):
            x[0] = x0
            for k in range(order):
                x[k + 1] = a * x[k] / (k + 1)
            x0 = t_horner(x, order, δt)
            output(x0, d0, d0, step * δt)


main()
