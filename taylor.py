#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from gmpy2 import mpfr, sqrt, exp, sinh, cosh, sin, cos, tanh, tan, log, asin, acos, atan

# noinspection PyArgumentList
to_mpfr = lambda x: mpfr(str(x)) if isinstance(x, (float, int)) else mpfr(x)

def t_jet(n, value=0):
    jet = [to_mpfr(0)] * n
    jet[0] = to_mpfr(value)
    return jet

def t_horner(jet, n, h):
    result = jet[n]
    for i in range(n - 1, -1, -1):
        result = result * h + jet[i]
    return result

t_prod = lambda u, v, k: sum(u[j] * v[k - j] for j in range(k + 1))

t_quot = lambda q, u, v, k: (u[k] - sum(v[j] * q[k - j] for j in range(1, k + 1))) / v[0]

ddot = lambda v, u, k: sum(j * u[j] * v[k - j] for j in range(1, k)) / k

t_sqrt = lambda r, u, k: sqrt(u[0]) if k == 0 else (u[k] / 2 - ddot(r, r, k)) / r[0]

t_exp = lambda e, u, k: exp(u[0]) if k == 0 else e[0] * u[k] + ddot(e, u, k)

def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    else:
        sn = c[0] * u[k] + ddot(c, u, k)
        cn = s[0] * u[k] + ddot(s, u, k)
        return sn, cn if hyp else - cn

def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1 - tanh(u[0])**2) if hyp else (tan(u[0]), tan(u[0])**2 + 1)
    else:
        tn = s2[0] * u[k] + ddot(s2, u, k)
        s2 = 2 * (t[0] * tn + ddot(t, t, k))
        return tn, - s2 if hyp else s2

t_pwr = lambda p, u, a, k: u[0]**a if k == 0 else (a * (p[0] * u[k] + ddot(p, u, k)) - ddot(u, p, k)) / u[0]

t_ln = lambda l, u, k: log(u[0]) if k == 0 else (u[k] - ddot(u, l, k)) / u[0]

t_atan = lambda h, v, u, k: (atan(u[0]), 1 + u[0]**2) if k == 0 else ((u[k] - ddot(v, h, k)) / v[0], 2 * (u[0] * u[k] + ddot(u, u, k)))

def _arc(h, v, u, k):
    h_jet = (u[k] - ddot(v, h, k)) / v[0]
    v_jet = u[0] * h_jet + ddot(u, h, k)
    return h_jet, - v_jet

t_asin = lambda h, v, u, k: (asin(u[0]), cos(asin(u[0]))) if k == 0 else _arc(h, v, u, k)

t_acos = lambda h, v, u, k: (acos(u[0]), - sin(acos(u[0]))) if k == 0 else _arc(h, v, u, k)

print(__name__ + " module loaded", file=stderr)
