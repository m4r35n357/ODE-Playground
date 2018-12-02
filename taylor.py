
from math import sin, cos, tan, exp, sinh, cosh, tanh, log, sqrt


def jet_0(n):
    return [0.0 for _ in range(n)]


def jet_c(constant, n):
    jet = jet_0(n)
    jet[0] = constant
    return jet


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
    assert v[0] != 0.0
    tq = 0.0
    for j in range(1, k + 1):
        tq += v[j] * q[k - j]
    return (u[k] - tq) / v[0]


def t_sqrt(r, u, k):
    assert u[0] > 0.0
    if k == 0:
        return sqrt(u[0])
    else:
        rt = 0.0
        if k % 2 == 1:
            for j in range(1, (k - 1) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt) / (2.0 * r[0])
        else:
            for j in range(1, (k - 2) // 2 + 1):
                rt += r[j] * r[k - j]
            return (u[k] - 2.0 * rt - r[k // 2]**2) / (2.0 * r[0])


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
    assert u[0] != 0.0
    if k == 0:
        return u[0]**a
    else:
        return (a * (p[0] * u[k] + ddot(p, u, k)) - ddot(u, p, k)) / u[0]


def t_ln(l, u, k):
    assert u[0] > 0.0
    if k == 0:
        return log(u[0])
    else:
        return (u[k] - ddot(u, l, k)) / u[0]
