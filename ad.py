#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv, stderr
from math import fsum, sin, cos, sinh, cosh, tan, tanh, exp, log, sqrt
from collections import namedtuple

class Components(namedtuple('ParametersType', ['x', 'y', 'z'])):
    pass

class Context:
    places = 3

def output(x, y, z, t):
    print(f'{x:+.{Context.places}e} {y:+.{Context.places}e} {z:+.{Context.places}e} {t:.5e}')

def t_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

def t_horner(jet, h):
    result = 0.0
    for term in reversed(jet):
        result = result * h + term
    return result

def t_abs(u, k):
    return - u[k] if u[0] < 0.0 else u[k]

def _cauchy(a, b, k, lower, upper):
    return fsum(a[j] * b[k - j] for j in range(lower, upper + 1))

def t_prod(u, v, k):
    return _cauchy(u, v, k, 0, k)

def t_quot(q, u, v, k):
    return u[0] / v[0] if k == 0 else (u[k] - _cauchy(q, v, k, 0, k - 1)) / v[0]

def t_inv(i, v, k):
    return 1.0 / v[0] if k == 0 else - _cauchy(i, v, k, 0, k - 1) / v[0]

def t_sqr(u, k):
    return _cauchy(u, u, k, 0, k)

def t_sqrt(r, u, k):
    return sqrt(u[0]) if k == 0 else 0.5 * (u[k] - _cauchy(r, r, k, 1, k - 1)) / r[0]

def _d_cauchy(h, u, k, lower, upper, factor=1.0):
    return factor * fsum(h[j] * (k - j) * u[k - j] for j in range(lower, upper + 1)) / k

def t_exp(e, u, k):
    return exp(u[0]) if k == 0 else _d_cauchy(e, u, k, 0, k - 1)

def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    s[k] = _d_cauchy(c, u, k, 0, k - 1)
    c[k] = _d_cauchy(s, u, k, 0, k - 1, -1.0 if not hyp else 1.0)
    return s[k], c[k]

def t_tan_sec2(t, s, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyp else (tan(u[0]), 1.0 + tan(u[0])**2)
    t[k] = _d_cauchy(s, u, k, 0, k - 1)
    s[k] = _d_cauchy(t, t, k, 0, k - 1, 2.0 if not hyp else -2.0)
    return t[k], s[k]

def t_pwr(p, u, a, k):
    return u[0]**a if k == 0 else (_d_cauchy(p, u, k, 0, k - 1, a) - _d_cauchy(u, p, k, 1, k - 1)) / u[0]

def t_ln(l, u, k):
    return log(u[0]) if k == 0 else (u[k] - _d_cauchy(u, l, k, 1, k - 1)) / u[0]

def tsm(ode, get_p, get_i):
    Context.places, n, δt, n_steps = argv[1], int(argv[3]), float(argv[4]), int(argv[5])  # controls
    x, y, z = t_jet(n + 1, float(argv[6])), t_jet(n + 1, float(argv[7])), t_jet(n + 1, float(argv[8]))  # initial values
    p = get_p(n)
    i = None if get_i is None else get_i(n)
    output(x[0], y[0], z[0], 0.0)
    for step in range(1, n_steps + 1):
        for k in range(n):
            c = ode(x, y, z, p, i, k)
            x[k + 1] = c.x / (k + 1)
            y[k + 1] = c.y / (k + 1)
            z[k + 1] = c.z / (k + 1)
        x[0], y[0], z[0] = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
        output(x[0], y[0], z[0], step * δt)

def rk4(ode, get_p):
    Context.places, interval, δt, n_steps = argv[1], int(argv[3]), float(argv[4]), int(argv[5])  # controls
    x, y, z = float(argv[6]), float(argv[7]), float(argv[8])  # initial values
    p = get_p()
    output(x, y, z, 0.0)
    for step in range(1, n_steps + 1):
        k1 = ode(x, y, z, p)
        k2 = ode(x + 0.5 * k1.x * δt, y + 0.5 * k1.y * δt, z + 0.5 * k1.z * δt, p)
        k3 = ode(x + 0.5 * k2.x * δt, y + 0.5 * k2.y * δt, z + 0.5 * k2.z * δt, p)
        k4 = ode(x + k3.x * δt, y + k3.y * δt, z + k3.z * δt, p)
        x += δt * (k1.x + 2.0 * (k2.x + k3.x) + k4.x) / 6.0
        y += δt * (k1.y + 2.0 * (k2.y + k3.y) + k4.y) / 6.0
        z += δt * (k1.z + 2.0 * (k2.z + k3.z) + k4.z) / 6.0
        if step % interval == 0:
            output(x, y, z, δt * step)

print(f'{__name__} module loaded', file=stderr)