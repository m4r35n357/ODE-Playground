#!/usr/bin/env python3
from collections import namedtuple
from math import sin, cos, tan, exp, sinh, cosh, tanh, sqrt
from sys import stderr


TrigTuple = namedtuple('TrigTupleType', ['pri', 'sec'])


def t_stepper(args):
    return int(args[1]), float(args[2]), int(args[3])


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


def t_sqrt(r, u, k):
    if k == 0:
        return sqrt(u[0])
    else:
        rt = 0.0
        if k % 2 == 1:
            for j in range(1, (k - 1) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt) / (2.0 * u[0])
        else:
            for j in range(1, (k - 2) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt - r[k // 2]**2) / (2.0 * u[0])


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


def t_sin_cos(s, c, u, k, hyperbolic=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyperbolic else (sin(u[0]), cos(u[0]))
    else:
        return c[0] * u[k] + ddot(c, u, k),\
               s[0] * u[k] + ddot(s, u, k) if hyperbolic else - (s[0] * u[k] - ddot(s, u, k))


def t_tan_sec2(t, s2, u, k, hyperbolic=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyperbolic else (tan(u[0]), tan(u[0])**2 + 1.0)
    else:
        t[k] = s2[0] * u[k] + ddot(s2, u, k)
        s2[k] = 2.0 * (t[0] * t[k] + ddot(t, t, k))
        return t[k], - s2[k] if hyperbolic else s2[k]


def t_pwr(p, u, a, k):
    if k == 0:
        return u[0]**a
    else:
        return (a * (p[0] * u[k] + ddot(p, u, k)) - ddot(u, p, k)) / u[0]


def t_ln(l, u, k):
    if k == 0:
        return exp(u[0])
    else:
        return (u[k] - ddot(u, l, k)) / u[0]


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
    return TrigTuple(pri=sin_jet, sec=cos_jet)


def ad_tan_sec2(u, hyperbolic=False):
    n = len(u)
    tan_jet = jet_0(n)
    sec2_jet = jet_0(n)
    for k in range(n):
        tan_jet[k], sec2_jet[k] = t_tan_sec2(tan_jet, sec2_jet, u, k, hyperbolic)
    return TrigTuple(pri=tan_jet, sec=sec2_jet)


def ad_ln(u):
    n = len(u)
    ln_jet = jet_0(n)
    for k in range(n):
        ln_jet[k] = t_ln(ln_jet, u, k)
    return ln_jet


def ad_sqrt(u):
    n = len(u)
    sqrt_jet = jet_0(n)
    for k in range(n):
        sqrt_jet[k] = t_sqrt(sqrt_jet, u, k)
    return sqrt_jet


def ad_pwr(u, a):
    n = len(u)
    pwr_jet = jet_0(n)
    for k in range(n):
        pwr_jet[k] = t_pwr(pwr_jet, u, a, k)
    return pwr_jet


print(__name__ + " module loaded", file=stderr)
