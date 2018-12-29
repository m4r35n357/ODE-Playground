#!/usr/bin/env python3

#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import argv
from gmpy2 import get_context, log, sqrt
get_context().precision = int(int(argv[2]) * log(10.0) / log(2.0))
from taylor import t_jet, t_horner, t_prod, t_sin_cos, t_tan_sec2, t_sqr, to_mpfr


def print_output(x, y, z, t):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, t))


def main():
    model = argv[1]
    n = int(argv[3])
    h = to_mpfr(argv[4])
    steps = int(argv[5])
    d0 = to_mpfr(0)

    x0, y0, z0 = to_mpfr(argv[6]), to_mpfr(argv[7]), to_mpfr(argv[8])
    x, y, z = t_jet(n + 1), t_jet(n + 1), t_jet(n + 1)

    if model == "lorenz":
        #  Example: ./tsm-mp.py lorenz 16 10 .01 3000 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        sigma, rho, beta = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]) / to_mpfr(argv[12])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = sigma * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (rho * x[k] - t_prod(x, z, k) - y[k]) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - beta * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "lu":
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "chen":
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = a * (y[k] - x[k]) / (k + 1)
                y[k + 1] = ((c - a) * x[k] + c * y[k] - t_prod(x, z, k)) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "rossler":
        a, b, c = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        b_ = t_jet(n, b)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = - (y[k] + z[k]) / (k + 1)
                y[k + 1] = (x[k] + a * y[k]) / (k + 1)
                z[k + 1] = (b_[k] + t_prod(x, z, k) - c * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "bouali":
        #  Example: ./tsm-mp.py bouali 80 40 0.02 50001 1 1 0 3 2.2 1 .01 | ./plotPi3d.py
        a, b, c, d = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]), to_mpfr(argv[12])
        w4 = t_jet(n)
        w5 = t_jet(n)
        jet1 = t_jet(n, 1)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                w4[k] = jet1[k] - y[k]
                w5[k] = jet1[k] - t_sqr(x, k)
                x[k + 1] = (a * t_prod(x, w4, k) - b * z[k]) / (k + 1)
                y[k + 1] = - c * t_prod(y, w5, k) / (k + 1)
                z[k + 1] = d * x[k] / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "thomas":
        #  Example: ./tsm-mp.py thomas 16 10 0.1 30000 1 0 0 .19 | ./plotPi3d.py
        b = to_mpfr(argv[9])
        sx, cx = t_jet(n), t_jet(n)
        sy, cy = t_jet(n), t_jet(n)
        sz, cz = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                sy[k], cy[k] = t_sin_cos(sy, cy, y, k)
                sz[k], cz[k] = t_sin_cos(sz, cz, z, k)
                x[k + 1] = (sy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (sz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (sx[k] - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "st":
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        sx, cx = t_jet(n), t_jet(n)
        sy, cy = t_jet(n), t_jet(n)
        sz, cz = t_jet(n), t_jet(n)
        ax, ay, az = t_jet(n), t_jet(n), t_jet(n)
        sax, cax = t_jet(n), t_jet(n)
        say, cay = t_jet(n), t_jet(n)
        saz, caz = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
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
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "rf":
        a, g = to_mpfr(argv[9]), to_mpfr(argv[10])
        w_a = t_jet(n)
        w_b = t_jet(n)
        w_c = t_jet(n)
        jet1 = t_jet(n, 1)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                w_x2_1 = t_sqr(x, k) - jet1[k]
                w_a[k] = z[k] + w_x2_1
                w_b[k] = 3 * z[k] - w_x2_1
                w_c[k] = a + t_prod(x, y, k)
                x[k + 1] = (t_prod(y, w_a, k) + g * x[k]) / (k + 1)
                y[k + 1] = (t_prod(x, w_b, k) + g * y[k]) / (k + 1)
                z[k + 1] = - 2.0 * t_prod(z, w_c, k) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "sprott":
        #  Example: ./tsm-mp.py sprott 16 10 0.1 30001 1 0 0 | ./plotPi3d.py
        w1 = t_jet(n, 1.0)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = (y[k] + 2 * t_prod(x, y, k) + t_prod(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - 2 * t_sqr(x, k) + t_prod(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - t_sqr(x, k) - t_sqr(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "sj":
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        w_b = t_jet(n, b)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - x[k] + t_prod(y, z, k) / (k + 1)
                z[k + 1] = (z[k] + a * t_sqr(x, k) - t_sqr(y, k) - w_b[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "halvorsen":
        #  Example: ./tsm-mp.py halvorsen 16 10 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        a = to_mpfr(argv[9])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = - (a * x[k] + 4 * y[k] + 4 * z[k] + t_sqr(y, k)) / (k + 1)
                y[k + 1] = - (a * y[k] + 4 * z[k] + 4 * x[k] + t_sqr(z, k)) / (k + 1)
                z[k + 1] = - (a * z[k] + 4 * x[k] + 4 * y[k] + t_sqr(x, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "nh":
        a_ = t_jet(n, to_mpfr(argv[9]))
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = (t_prod(y, z, k) - x[k]) / (k + 1)
                z[k + 1] = (a_[k] - t_sqr(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "rucklidge":
        a, b = to_mpfr(argv[9]), to_mpfr(argv[10])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = (a * y[k] - b * x[k] - t_prod(y, z, k)) / (k + 1)
                y[k + 1] = x[k] / (k + 1)
                z[k + 1] = (t_sqr(y, k) - z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "damped":
        c1, c2 = to_mpfr(argv[9]), to_mpfr(argv[10])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - (c1 * x[k] + c2 * y[k]) / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "pendulum":
        w = sqrt(to_mpfr(argv[9]) / to_mpfr(argv[10]))
        sx, cx = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - w * sx[k] / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "volterra":
        a, b, c, d = to_mpfr(argv[8]), to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                wxy = t_prod(x, y, k)
                x[k + 1] = (a * x[k] - c * wxy) / (k + 1)
                y[k + 1] = (d * wxy - b * y[k]) / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)


main()
