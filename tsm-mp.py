#!/usr/bin/env python3

#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from math import sqrt
from sys import stderr, argv
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

    x, y, z = float(argv[6]), float(argv[7]), float(argv[8])
    cx, cy, cz = t_jet(n + 1), t_jet(n + 1), t_jet(n + 1)

    if model == "lorenz":
        #  Example: ./tsm-mp.py lorenz 16 10 .01 100001 -15.8 -17.48 35.64 10 28 8 3 | ./plotPi3d.py
        s, r, b = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = s * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = (r * cx[k] - t_prod(cx, cz, k) - cy[k]) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "lu":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = a * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = (c * cy[k] - t_prod(cx, cz, k)) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "chen":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = a * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = ((c - a) * cx[k] + c * cy[k] - t_prod(cx, cz, k)) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "rossler":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        b_ = t_jet(n, b)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = - (cy[k] + cz[k]) / (k + 1)
                cy[k + 1] = (cx[k] + a * cy[k]) / (k + 1)
                cz[k + 1] = (b_[k] + t_prod(cx, cz, k) - c * cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "bouali":
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        w4 = t_jet(n)
        w5 = t_jet(n)
        jet1 = t_jet(n, d1)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                w4[k] = jet1[k] - cy[k]
                w5[k] = jet1[k] - t_sqr(cx, k)
                cx[k + 1] = (a * t_prod(cx, w4, k) - b * cz[k]) / (k + 1)
                cy[k + 1] = - c * t_prod(cy, w5, k) / (k + 1)
                cz[k + 1] = d * cx[k] / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "thomas":
        #  Example: ./tsm-mp.py thomas 16 10 0.1 30001 1 0 0 .19 | ./plotPi3d.py
        b = float(argv[9])
        wsx, wcx = t_jet(n), t_jet(n)
        wsy, wcy = t_jet(n), t_jet(n)
        wsz, wcz = t_jet(n), t_jet(n)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, cx, k)
                wsy[k], wcy[k] = t_sin_cos(wsy, wcy, cy, k)
                wsz[k], wcz[k] = t_sin_cos(wsz, wcz, cz, k)
                cx[k + 1] = (wsy[k] - b * cx[k]) / (k + 1)
                cy[k + 1] = (wsz[k] - b * cy[k]) / (k + 1)
                cz[k + 1] = (wsx[k] - b * cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "st":
        a, b = float(argv[9]), float(argv[10])
        wsx, wcx = t_jet(n), t_jet(n)
        wsy, wcy = t_jet(n), t_jet(n)
        wsz, wcz = t_jet(n), t_jet(n)
        wax, way, waz = t_jet(n), t_jet(n), t_jet(n)
        wsax, wcax = t_jet(n), t_jet(n)
        wsay, wcay = t_jet(n), t_jet(n)
        wsaz, wcaz = t_jet(n), t_jet(n)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                wsx[k], wcx[k] = t_tan_sec2(wsx, wcx, cx, k)
                wsy[k], wcy[k] = t_tan_sec2(wsy, wcy, cy, k)
                wsz[k], wcz[k] = t_tan_sec2(wsz, wcz, cz, k)
                wax[k], way[k], waz[k] = a * cx[k], a * cy[k], a * cz[k]
                wsax[k], wcax[k] = t_sin_cos(wsax, wcax, wax, k)
                wsay[k], wcay[k] = t_sin_cos(wsay, wcay, way, k)
                wsaz[k], wcaz[k] = t_sin_cos(wsaz, wcaz, waz, k)
                cx[k + 1] = (wsay[k] - b * wsx[k]) / (k + 1)
                cy[k + 1] = (wsaz[k] - b * wsy[k]) / (k + 1)
                cz[k + 1] = (wsax[k] - b * wsz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "rf":
        a, g = float(argv[9]), float(argv[10])
        w_a = t_jet(n)
        w_b = t_jet(n)
        w_c = t_jet(n)
        jet1 = t_jet(n, d1)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                w_x2_1 = t_sqr(cx, k) - jet1[k]
                w_a[k] = cz[k] + w_x2_1
                w_b[k] = d3 * cz[k] - w_x2_1
                w_c[k] = a + t_prod(cx, cy, k)
                cx[k + 1] = (t_prod(cy, w_a, k) + g * cx[k]) / (k + 1)
                cy[k + 1] = (t_prod(cx, w_b, k) + g * cy[k]) / (k + 1)
                cz[k + 1] = - 2.0 * t_prod(cz, w_c, k) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "sprott":
        #  Example: ./tsm-mp.py sprott 16 10 0.1 30001 1 0 0 | ./plotPi3d.py
        w1 = t_jet(n, 1.0)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = (cy[k] + d2 * t_prod(cx, cy, k) + t_prod(cx, cz, k)) / (k + 1)
                cy[k + 1] = (w1[k] - d2 * t_sqr(cx, k) + t_prod(cy, cz, k)) / (k + 1)
                cz[k + 1] = (cx[k] - t_sqr(cx, k) - t_sqr(cy, k)) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "sj":
        a, b = float(argv[9]), float(argv[10])
        w_b = t_jet(n, b)
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - cx[k] + t_prod(cy, cz, k) / (k + 1)
                cz[k + 1] = (cz[k] + a * t_sqr(cx, k) - t_sqr(cy, k) - w_b[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "halvorsen":
        #  Example: ./tsm-mp.py halvorsen 16 10 .01 100001 1 0 0 1.4 | ./plotPi3d.py
        a = float(argv[9])
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = - (a * cx[k] + d4 * cy[k] + d4 * cz[k] + t_sqr(cy, k)) / (k + 1)
                cy[k + 1] = - (a * cy[k] + d4 * cz[k] + d4 * cx[k] + t_sqr(cz, k)) / (k + 1)
                cz[k + 1] = - (a * cz[k] + d4 * cx[k] + d4 * cy[k] + t_sqr(cx, k)) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "nh":
        a_ = t_jet(n, float(argv[9]))
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = (t_prod(cy, cz, k) - cx[k]) / (k + 1)
                cz[k + 1] = (a_[k] - t_sqr(cy, k)) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "rucklidge":
        a, b = float(argv[9]), float(argv[10])
        for step in range(1, steps):
            print_output(x, y, z, step * h)
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = (a * cy[k] - b * cx[k] - t_prod(cy, cz, k)) / (k + 1)
                cy[k + 1] = cx[k] / (k + 1)
                cz[k + 1] = (t_sqr(cy, k) - cz[k]) / (k + 1)
            x, y, z = t_horner(cx, n, h), t_horner(cy, n, h), t_horner(cz, n, h)
    elif model == "damped":
        c1, c2 = float(argv[9]), float(argv[10])
        for step in range(1, steps):
            print_output(x, y, d0, step * h)
            cx[0], cy[0] = x, y
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - (c1 * cx[k] + c2 * cy[k]) / (k + 1)
            x, y = t_horner(cx, n, h), t_horner(cy, n, h)
    elif model == "pendulum":
        w = sqrt(float(argv[9]) / float(argv[10]))
        wsx, wcx = t_jet(n), t_jet(n)
        for step in range(1, steps):
            print_output(x, y, d0, step * h)
            cx[0], cy[0] = x, y
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, cx, k)
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - w * wsx[k] / (k + 1)
            x, y = t_horner(cx, n, h), t_horner(cy, n, h)
    elif model == "volterra":
        a, b, c, d = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print_output(x, y, d0, step * h)
            cx[0], cy[0] = x, y
            for k in range(n):
                wxy = t_prod(cx, cy, k)
                cx[k + 1] = (a * cx[k] - c * wxy) / (k + 1)
                cy[k + 1] = (d * wxy - b * cy[k]) / (k + 1)
            x, y = t_horner(cx, n, h), t_horner(cy, n, h)


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
