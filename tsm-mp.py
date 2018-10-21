#!/usr/bin/env python3

from math import sin, cos, tan, sqrt
from sys import stderr, argv

D0 = float('0.0')
D05 = float('0.5')
D1 = float('1.0')
D2 = float('2.0')
D3 = float('3.0')
D4 = float('4.0')


def jet_0(n):
    return [D0 for _ in range(n)]


def jet_c(n, constant):
    jet = jet_0(n)
    jet[0] = constant
    return jet


def t_prod(v, u, k):
    tp = D0
    for j in range(k + 1):
        tp += u[j] * v[k - j]
    return tp


def t_chain(df_du, u, k):
    df_dx = D0
    for j in range(k):
        df_dx += df_du[j] * (k - j) * u[k - j]
    return df_dx / k


def t_sin_cos(s, c, u, k):
    if k == 0:
        return sin(u[0]), cos(u[0])
    else:
        return t_chain(c, u, k), - t_chain(s, u, k)


def t_tan_sec2(t, s2, u, k):
    if k == 0:
        return tan(u[0]), tan(u[0])**2 + D1
    else:
        return t_chain(s2, u, k), D2 * t_chain(t, t, k)


def horner(jet, n, h):
    result = jet[n]
    for i in range(n - 1, -1, -1):
        result = result * h + jet[i]
    return result


def main():
    model = argv[1]
    n = int(argv[3])
    h = float(argv[4])
    steps = int(argv[5])

    x = float(argv[6])
    y = float(argv[7])
    z = float(argv[8])

    cx = jet_0(n + 1)
    cy = jet_0(n + 1)
    cz = jet_0(n + 1)
    if model == "lorenz":
        s, r, b = float(argv[9]), float(argv[10]), float(argv[11]) / float(argv[12])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = s * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = (r * cx[k] - t_prod(cx, cz, k) - cy[k]) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "lu":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = a * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = (c * cy[k] - t_prod(cx, cz, k)) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "chen":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = a * (cy[k] - cx[k]) / (k + 1)
                cy[k + 1] = ((c - a) * cx[k] + c * cy[k] - t_prod(cx, cz, k)) / (k + 1)
                cz[k + 1] = (t_prod(cx, cy, k) - b * cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "rossler":
        a, b, c = float(argv[9]), float(argv[10]), float(argv[11])
        b_ = jet_c(n, b)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = - (cy[k] + cz[k]) / (k + 1)
                cy[k + 1] = (cx[k] + a * cy[k]) / (k + 1)
                cz[k + 1] = (b_[k] + t_prod(cx, cz, k) - c * cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "bouali":
        a, b, c, d = float(argv[9]), float(argv[10]), float(argv[11]), float(argv[12])
        w4 = jet_0(n)
        w5 = jet_0(n)
        jet1 = jet_c(n, D1)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                w4[k] = jet1[k] - cy[k]
                w5[k] = jet1[k] - t_prod(cx, cx, k)
                cx[k + 1] = (a * t_prod(cx, w4, k) - b * cz[k]) / (k + 1)
                cy[k + 1] = - c * t_prod(cy, w5, k) / (k + 1)
                cz[k + 1] = d * cx[k] / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "thomas":
        b = float(argv[9])
        wsx, wcx = jet_0(n), jet_0(n)
        wsy, wcy = jet_0(n), jet_0(n)
        wsz, wcz = jet_0(n), jet_0(n)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, cx, k)
                wsy[k], wcy[k] = t_sin_cos(wsy, wcy, cy, k)
                wsz[k], wcz[k] = t_sin_cos(wsz, wcz, cz, k)
                cx[k + 1] = (wsy[k] - b * cx[k]) / (k + 1)
                cy[k + 1] = (wsz[k] - b * cy[k]) / (k + 1)
                cz[k + 1] = (wsx[k] - b * cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "st":
        a, b = float(argv[9]), float(argv[10])
        wsx, wcx = jet_0(n), jet_0(n)
        wsy, wcy = jet_0(n), jet_0(n)
        wsz, wcz = jet_0(n), jet_0(n)
        wax, way, waz = jet_0(n), jet_0(n), jet_0(n)
        wsax, wcax = jet_0(n), jet_0(n)
        wsay, wcay = jet_0(n), jet_0(n)
        wsaz, wcaz = jet_0(n), jet_0(n)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
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
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "rf":
        a, g = float(argv[9]), float(argv[10])
        w_a = jet_0(n)
        w_b = jet_0(n)
        w_c = jet_0(n)
        jet1 = jet_c(n, D1)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                w_x2_1 = t_prod(cx, cx, k) - jet1[k]
                w_a[k] = cz[k] + w_x2_1
                w_b[k] = D3 * cz[k] - w_x2_1
                w_c[k] = a + t_prod(cx, cy, k)
                cx[k + 1] = (t_prod(cy, w_a, k) + g * cx[k]) / (k + 1)
                cy[k + 1] = (t_prod(cx, w_b, k) + g * cy[k]) / (k + 1)
                cz[k + 1] = - 2.0 * t_prod(cz, w_c, k) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "sj":
        a, b = float(argv[9]), float(argv[10])
        w_b = jet_c(n, b)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - cx[k] + t_prod(cy, cz, k) / (k + 1)
                cz[k + 1] = (cz[k] + a * t_prod(cx, cx, k) - t_prod(cy, cy, k) - w_b[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "halvorsen":
        a = float(argv[9])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = - (a * cx[k] + D4 * cy[k] + D4 * cz[k] + t_prod(cy, cy, k)) / (k + 1)
                cy[k + 1] = - (a * cy[k] + D4 * cz[k] + D4 * cx[k] + t_prod(cz, cz, k)) / (k + 1)
                cz[k + 1] = - (a * cz[k] + D4 * cx[k] + D4 * cy[k] + t_prod(cx, cx, k)) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "nh":
        a_ = jet_c(n, float(argv[9]))
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = (t_prod(cy, cz, k) - cx[k]) / (k + 1)
                cz[k + 1] = (a_[k] - t_prod(cy, cy, k)) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "rucklidge":
        a, b = float(argv[9]), float(argv[10])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, z, step * h))
            cx[0], cy[0], cz[0] = x, y, z
            for k in range(n):
                cx[k + 1] = (a * cy[k] - b * cx[k] - t_prod(cy, cz, k)) / (k + 1)
                cy[k + 1] = cx[k] / (k + 1)
                cz[k + 1] = (t_prod(cy, cy, k) - cz[k]) / (k + 1)
            x, y, z = horner(cx, n, h), horner(cy, n, h), horner(cz, n, h)
    elif model == "damped":
        c1, c2 = float(argv[9]), float(argv[10])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, D0, step * h))
            cx[0], cy[0] = x, y
            for k in range(n):
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - (c1 * cx[k] + c2 * cy[k]) / (k + 1)
            x, y = horner(cx, n, h), horner(cy, n, h)
    elif model == "pendulum":
        w = sqrt(float(argv[9]) / float(argv[10]))
        wsx, wcx = jet_0(n), jet_0(n)
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, D0, step * h))
            cx[0], cy[0] = x, y
            for k in range(n):
                wsx[k], wcx[k] = t_sin_cos(wsx, wcx, cx, k)
                cx[k + 1] = cy[k] / (k + 1)
                cy[k + 1] = - w * wsx[k] / (k + 1)
            x, y = horner(cx, n, h), horner(cy, n, h)
    elif model == "volterra":
        a, b, c, d = float(argv[8]), float(argv[9]), float(argv[10]), float(argv[11])
        for step in range(1, steps):
            print("{:.9e} {:.9e} {:.9e} {:.5e}".format(x, y, D0, step * h))
            cx[0], cy[0] = x, y
            for k in range(n):
                wxy = t_prod(cx, cy, k)
                cx[k + 1] = (a * cx[k] - c * wxy) / (k + 1)
                cy[k + 1] = (d * wxy - b * cy[k]) / (k + 1)
            x, y = horner(cx, n, h), horner(cy, n, h)


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
