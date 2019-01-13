#!/usr/bin/env python3

#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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
        #  Example: ./tsm-mp.py rossler 16 10 0.01 150000 0.0 -6.78 0.02 .2 .2 5.7 | ./plotPi3d.py
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
        jet1, w4, w5 = t_jet(n, 1), t_jet(n), t_jet(n)
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
        #  Example: ./tsm-mp.py rf 16 10 .01 100001 .1 .1 .1 .2876 .1 | ./plotPi3d.py
        a, g = to_mpfr(argv[9]), to_mpfr(argv[10])
        jet1, w_a, w_b, w_c = t_jet(n, 1), t_jet(n), t_jet(n), t_jet(n)
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
                z[k + 1] = - 2 * t_prod(z, w_c, k) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "sprott":
        #  Example: ./tsm-mp.py sprott 16 10 0.1 30001 1 0 0 | ./plotPi3d.py
        w1 = t_jet(n, 1)
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
        #  Example: ./tsm-mp.py rucklidge 16 10 0.01 10001 1 0 0 6.7 2 | ./plotPi3d.py
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
        print_output(x0, y0, d0, d0)
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
        print_output(x0, y0, d0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                sx[k], cx[k] = t_sin_cos(sx, cx, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - w * sx[k] / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "volterra":
        #  Example: ./tsm-mp.py volterra 16 10 .01 2001 10 10 0 1 .5 .05 .02 | ./plotAnimated.py 1 0 80
        a, b, c, d = to_mpfr(argv[9]), to_mpfr(argv[10]), to_mpfr(argv[11]), to_mpfr(argv[12])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                wxy = t_prod(x, y, k)
                x[k + 1] = (a * x[k] - c * wxy) / (k + 1)
                y[k + 1] = (d * wxy - b * y[k]) / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "logistic":
        #  Example: ./tsm-mp.py logistic 16 10 0.1 10001 .6 0 0 .1 | ./plotXY.py 1 3 0
        a = to_mpfr(argv[9])
        w1, wa, wb = t_jet(n, 1), t_jet(n), t_jet(n)
        print_output(x0, d0, d0, d0)
        for step in range(1, steps + 1):
            x[0] = x0
            for k in range(n):
                wa[k] = a * x[k]
                wb[k] = w1[k] - x[k]
                x[k + 1] = t_prod(wa, wb, k) / (k + 1)
            x0 = t_horner(x, n, h)
            print_output(x0, d0, d0, step * h)
    elif model == "constant":
        #  Example: ./tsm-mp.py constant 16 10 0.1 10001 10 0 0 -.05 | ./plotXY.py 1 3 0
        a = to_mpfr(argv[9])
        print_output(x0, d0, d0, d0)
        for step in range(1, steps + 1):
            x[0] = x0
            for k in range(n):
                x[k + 1] = a * x[k] / (k + 1)
            x0 = t_horner(x, n, h)
            print_output(x0, d0, d0, step * h)


main()
