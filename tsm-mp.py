#!/usr/bin/env python3

#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from math import sqrt
from sys import argv
from taylor import t_jet, t_horner, t_prod, t_sin_cos, t_tan_sec2, t_sqr


def print_output(x, y, z, t):
    print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, t))


def main():
    model = argv[1]
    n = int(argv[3])
    h = float(argv[4])
    steps = int(argv[5])
    d0 = float('0.0')
    d1 = float('1.0')
    d2 = float('2.0')
    d3 = float('3.0')
    d4 = float('4.0')

    x0, y0, z0 = float(argv[6]), float(argv[7]), float(argv[8])
    x, y, z = t_jet(n + 1), t_jet(n + 1), t_jet(n + 1)

    if model == "lorenz":
        #  Example: ./tsm-mp.py lorenz 16 10 .01 100001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        s, r, b = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = s * (y[k] - x[k]) / (k + 1)
                y[k + 1] = (r * x[k] - t_prod(x, z, k) - y[k]) / (k + 1)
                z[k + 1] = (t_prod(x, y, k) - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "lu":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
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
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        w4 = t_jet(n)
        w5 = t_jet(n)
        jet1 = t_jet(n, d1)
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
        #  Example: ./tsm-mp.py thomas 16 10 0.1 30001 1 0 0 .19 | ./plotPi3d.py
        b = float(argv[9])
        wsx, wcx = t_jet(n), t_jet(n)
        wsy, wcy = t_jet(n), t_jet(n)
        wsz, wcz = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, x, k)
                wsy[k], wcy[k] = t_sin_cos(wsy, wcy, y, k)
                wsz[k], wcz[k] = t_sin_cos(wsz, wcz, z, k)
                x[k + 1] = (wsy[k] - b * x[k]) / (k + 1)
                y[k + 1] = (wsz[k] - b * y[k]) / (k + 1)
                z[k + 1] = (wsx[k] - b * z[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "st":
        a, b = float(argv[9]), float(argv[10])
        wsx, wcx = t_jet(n), t_jet(n)
        wsy, wcy = t_jet(n), t_jet(n)
        wsz, wcz = t_jet(n), t_jet(n)
        wax, way, waz = t_jet(n), t_jet(n), t_jet(n)
        wsax, wcax = t_jet(n), t_jet(n)
        wsay, wcay = t_jet(n), t_jet(n)
        wsaz, wcaz = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                wsx[k], wcx[k] = t_tan_sec2(wsx, wcx, x, k)
                wsy[k], wcy[k] = t_tan_sec2(wsy, wcy, y, k)
                wsz[k], wcz[k] = t_tan_sec2(wsz, wcz, z, k)
                wax[k], way[k], waz[k] = a * x[k], a * y[k], a * z[k]
                wsax[k], wcax[k] = t_sin_cos(wsax, wcax, wax, k)
                wsay[k], wcay[k] = t_sin_cos(wsay, wcay, way, k)
                wsaz[k], wcaz[k] = t_sin_cos(wsaz, wcaz, waz, k)
                x[k + 1] = (wsay[k] - b * wsx[k]) / (k + 1)
                y[k + 1] = (wsaz[k] - b * wsy[k]) / (k + 1)
                z[k + 1] = (wsax[k] - b * wsz[k]) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "rf":
        a, g = float(argv[9]), float(argv[10])
        w_a = t_jet(n)
        w_b = t_jet(n)
        w_c = t_jet(n)
        jet1 = t_jet(n, d1)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                w_x2_1 = t_sqr(x, k) - jet1[k]
                w_a[k] = z[k] + w_x2_1
                w_b[k] = d3 * z[k] - w_x2_1
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
                x[k + 1] = (y[k] + d2 * t_prod(x, y, k) + t_prod(x, z, k)) / (k + 1)
                y[k + 1] = (w1[k] - d2 * t_sqr(x, k) + t_prod(y, z, k)) / (k + 1)
                z[k + 1] = (x[k] - t_sqr(x, k) - t_sqr(y, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "sj":
        a, b = float(argv[9]), float(argv[10])
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
        a = float(argv[9])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0], z[0] = x0, y0, z0
            for k in range(n):
                x[k + 1] = - (a * x[k] + d4 * y[k] + d4 * z[k] + t_sqr(y, k)) / (k + 1)
                y[k + 1] = - (a * y[k] + d4 * z[k] + d4 * x[k] + t_sqr(z, k)) / (k + 1)
                z[k + 1] = - (a * z[k] + d4 * x[k] + d4 * y[k] + t_sqr(x, k)) / (k + 1)
            x0, y0, z0 = t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h)
            print_output(x0, y0, z0, step * h)
    elif model == "nh":
        a_ = t_jet(n, float(argv[9]))
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
        a, b = float(argv[9]), float(argv[10])
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
        c1, c2 = float(argv[9]), float(argv[10])
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - (c1 * x[k] + c2 * y[k]) / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "pendulum":
        w = sqrt(float(argv[9]) / float(argv[10]))
        wsx, wcx = t_jet(n), t_jet(n)
        print_output(x0, y0, z0, d0)
        for step in range(1, steps + 1):
            x[0], y[0] = x0, y0
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, x, k)
                x[k + 1] = y[k] / (k + 1)
                y[k + 1] = - w * wsx[k] / (k + 1)
            x0, y0 = t_horner(x, n, h), t_horner(y, n, h)
            print_output(x0, y0, d0, step * h)
    elif model == "volterra":
        a, b, c, d = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
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
