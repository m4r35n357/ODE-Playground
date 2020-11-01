#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import argv, stderr
from math import fsum, sin, cos, sinh, cosh, tan, tanh, exp
from collections import namedtuple

class Components(namedtuple('ParametersType', ['x', 'y', 'z'])):
    pass

class Context:
    places = 3

def output(x, y, z, t, x_label="_", y_label="_", z_label="_"):
    print(f'{x:+.{Context.places}e} {y:+.{Context.places}e} {z:+.{Context.places}e} {t:.5e} {x_label} {y_label} {z_label}')

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

def t_sqr(u, k):
    return _cauchy(u, u, k, 0, k)

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

def tsm(ode, get_p, get_i):
    Context.places, n, δt, n_steps = argv[1], int(argv[2]), float(argv[3]), int(argv[4])  # controls
    x, y, z = t_jet(n + 1, float(argv[5])), t_jet(n + 1, float(argv[6])), t_jet(n + 1, float(argv[7]))  # initial values
    p = get_p(n)
    i = None if get_i is None else get_i(n)
    output(x[0], y[0], z[0], 0.0)
    for step in range(1, n_steps + 1):
        xdot, ydot, zdot = x[1], y[1], z[1]
        for k in range(n):
            c = ode(x, y, z, p, i, k)
            x[k + 1] = c.x / (k + 1)
            y[k + 1] = c.y / (k + 1)
            z[k + 1] = c.z / (k + 1)
        x[0], y[0], z[0] = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
        output(x[0], y[0], z[0], step * δt,
               ("x" if x[2] > 0.0 else "X") if x[1] * xdot < 0.0 else "_",
               ("y" if y[2] > 0.0 else "Y") if y[1] * ydot < 0.0 else "_",
               ("z" if z[2] > 0.0 else "Z") if z[1] * zdot < 0.0 else "_")

print(f'{__name__} module loaded', file=stderr)