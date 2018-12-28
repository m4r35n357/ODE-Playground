#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from gmpy2 import mpfr, sqrt, exp, sinh, cosh, sin, cos, tanh, tan, log, asin, acos, atan

# noinspection PyArgumentList
D0 = mpfr("0.0")
# noinspection PyArgumentList
D1 = mpfr("1.0")
# noinspection PyArgumentList
D2 = mpfr("2.0")


def t_jet(n, value=D0):
    jet = [D0 for _ in range(n)]
    jet[0] = value
    return jet


def t_horner(jet, n, h):
    result = jet[n]
    for i in range(n - 1, -1, -1):
        result = result * h + jet[i]
    return result


def t_prod(u, v, k):
    return sum(u[j] * v[k - j] for j in range(k + 1))


def t_quot(q, u, v, k):
    return (u[k] - sum(v[j] * q[k - j] for j in range(1, k + 1))) / v[0]


def _ddot(v, u, k):
    return sum(j * u[j] * v[k - j] for j in range(1, k)) / k


def t_sqr(u, k):
    if k == 0:
        return u[0]**2
    else:
        return D2 * (u[0] * u[k] + _ddot(u, u, k))


def t_sqrt(r, u, k):
    if k == 0:
        return sqrt(u[0])
    else:
        return (u[k] / D2 - _ddot(r, r, k)) / r[0]


def t_exp(e, u, k):
    if k == 0:
        return exp(u[0])
    else:
        return e[0] * u[k] + _ddot(e, u, k)


def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    else:
        sn = c[0] * u[k] + _ddot(c, u, k)
        cn = s[0] * u[k] + _ddot(s, u, k)
        return sn, cn if hyp else - cn


def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), D1 - tanh(u[0])**2) if hyp else (tan(u[0]), tan(u[0])**2 + D1)
    else:
        tn = s2[0] * u[k] + _ddot(s2, u, k)
        s2 = D2 * (t[0] * tn + _ddot(t, t, k))
        return tn, - s2 if hyp else s2


def t_pwr(p, u, a, k):
    if k == 0:
        return u[0]**a
    else:
        return (a * (p[0] * u[k] + _ddot(p, u, k)) - _ddot(u, p, k)) / u[0]


def t_ln(l, u, k):
    if k == 0:
        return log(u[0])
    else:
        return (u[k] - _ddot(u, l, k)) / u[0]


def t_atan(h, v, u, k):
    if k == 0:
        return atan(u[0]), D1 + u[0]**2
    else:
        return (u[k] - _ddot(v, h, k)) / v[0], D2 * (u[0] * u[k] + _ddot(u, u, k))


def _arc(h, v, u, k):
    h_jet = (u[k] - _ddot(v, h, k)) / v[0]
    v_jet = u[0] * h_jet + _ddot(u, h, k)
    return h_jet, - v_jet


def t_asin(h, v, u, k):
    if k == 0:
        h = asin(u[0])
        return h, cos(h)
    else:
        return _arc(h, v, u, k)


def t_acos(h, v, u, k):
    if k == 0:
        h = acos(u[0])
        return h, - sin(h)
    else:
        return _arc(h, v, u, k)


print(__name__ + " module loaded", file=stderr)
