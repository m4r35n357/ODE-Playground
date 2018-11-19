#!/usr/bin/env python3

from math import sin, cos, tan, exp, sinh, cosh, tanh
from sys import stderr


def jet_0(n):
    return [0.0 for _ in range(n)]


def jet_c(constant, n, diff=False):
    jet = jet_0(n)
    jet[0] = constant
    if diff:
        jet[1] = 1.0
    return jet


def derivatives(jet):
    n = len(jet)
    fac = 1
    for i in range(1, n):
        fac *= i
        jet[i] *= fac


def jet_output(jet):
    n = len(jet)
    print("{:14.6e} ".format(jet[0]), end='')
    for i in range(1, n):
        print("{:14.6e}".format(jet[i]), end='')
    print("")


def t_horner(jet, n, h):
    result = jet[n]
    for i in range(n - 1, -1, -1):
        result = result * h + jet[i]
    return result


def t_sqr(u, k):
    if k == 0:
        return u[0]**2
    else:
        sq = 0.0
        if k % 2 == 1:
            for j in range((k - 1) // 2 + 1):
                sq += u[j] * u[k - j]
            return 2.0 * sq
        else:
            for j in range((k - 2) // 2 + 1):
                sq += u[j] * u[k - j]
            return 2.0 * sq + u[k // 2]**2


def t_prod(u, v, k):
    tp = 0.0
    for j in range(k + 1):
        tp += u[j] * v[k - j]
    return tp


def t_quot(q, u, v, k):
    tq = 0.0
    for j in range(1, k + 1):
        tq += v[j] * q[k - j]
    return (u[k] - tq) / v[0]


def ddot(v, u, k):
    dd = 0.0
    for j in range(1, k):
        dd += j * u[j] * v[k - j]
    return dd / k


def t_exp(e, u, k):
    if k == 0:
        return exp(u[0])
    else:
        return e[0] * u[k] + ddot(e, u, k)


def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    else:
        sn = c[0] * u[k] + ddot(c, u, k)
        cn = s[0] * u[k] + ddot(s, u, k)
        return sn, cn if hyp else - cn


def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyp else (tan(u[0]), tan(u[0])**2 + 1.0)
    else:
        tn = s2[0] * u[k] + ddot(s2, u, k)
        sc2 = 2.0 * (t[0] * tn + ddot(t, t, k))
        return tn, - sc2 if hyp else sc2


def t_pwr(p, u, a, k):
    if k == 0:
        return u[0]**a
    else:
        return (a * (p[0] * u[k] + ddot(p, u, k)) - ddot(u, p, k)) / u[0]


def ad_scale(u, a):
    n = len(u)
    jet = jet_0(n)
    for k in range(n):
        jet[k] = a * u[k]
    return jet


def ad_plus(u, v):
    assert len(u) == len(v)
    n = len(u)
    jet = jet_0(n)
    for k in range(n):
        jet[k] = u[k] + v[k]
    return jet


def ad_minus(u, v):
    assert len(u) == len(v)
    n = len(u)
    jet = jet_0(n)
    for k in range(n):
        jet[k] = u[k] - v[k]
    return jet


def ad_sqr(u):
    n = len(u)
    sqr_jet = jet_0(n)
    for k in range(n):
        sqr_jet[k] = t_sqr(u, k)
    return sqr_jet


def ad_prod(u, v):
    assert len(u) == len(v)
    n = len(u)
    jet = jet_0(n)
    for k in range(n):
        jet[k] = t_prod(u, v, k)
    return jet


def ad_quot_inv(u, v):
    assert len(u) == len(v)
    n = len(u)
    q = jet_0(n)
    for k in range(n):
        q[k] = t_quot(q, u, v, k)
    return q


def ad_exp(u):
    n = len(u)
    exp_jet = jet_0(n)
    for k in range(n):
        exp_jet[k] = t_exp(exp_jet, u, k)
    return exp_jet


def ad_sin_cos(u, hyperbolic=False):
    n = len(u)
    sin_jet = jet_0(n)
    cos_jet = jet_0(n)
    for k in range(n):
        sin_jet[k], cos_jet[k] = t_sin_cos(sin_jet, cos_jet, u, k, hyperbolic)
    return sin_jet, cos_jet


def ad_tan_sec2(u, hyperbolic=False):
    n = len(u)
    tan_jet = jet_0(n)
    sec2_jet = jet_0(n)
    for k in range(n):
        tan_jet[k], sec2_jet[k] = t_tan_sec2(tan_jet, sec2_jet, u, k, hyperbolic)
    return tan_jet, sec2_jet


def ad_pwr(u, a):
    n = len(u)
    pwr_jet = jet_0(n)
    for k in range(n):
        pwr_jet[k] = t_pwr(pwr_jet, u, a, k)
    return pwr_jet


def ad_newton(model, x, target=0.0, tol=1.0e-12, max_it=100):
    f = [1.0, 0.0]
    delta = 1.0
    counter = 1
    while abs(f[0]) > tol or abs(delta) > tol:
        f = model(x, target)
        delta = - f[0] / f[1]
        x[0] += delta
        print("{:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, x[0], f[0] + target, delta))
        counter += 1
        if counter > max_it:
            break


if __name__ == "__main__":
    def fun(x, value):
        body = ad_sqr(x)
        return ad_minus(body, jet_c(value, len(x)))

    ad_newton(fun, jet_c(2.0, 2, diff=True), target=64.0)
    print("")

    from math import pi
    sine, cosine = ad_sin_cos(jet_c(pi / 3.0, 7, diff=True))
    derivatives(sine)
    jet_output(sine)
    derivatives(cosine)
    jet_output(cosine)
    print("")
    tangent, secant2 = ad_tan_sec2(jet_c(pi / 4.0, 7, diff=True))
    derivatives(tangent)
    jet_output(tangent)
    derivatives(secant2)
    jet_output(secant2)

else:
    print(__name__ + " module loaded", file=stderr)
