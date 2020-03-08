import system/ansi_c
import math

type
  Geometry* = enum
    trig, hyp

proc t_output* (x, y, z, t: float) =
    c_printf("%+.12e %+.12e %+.12e %+.6e\n", x, y, z, t);  #echo(x, ' ', y, ' ', z, ' ', t)

proc t_jet* (n: int): auto =
    return newSeq[float](n + 1)

proc t_jet_c* (n: int, value:float): auto =
    result = t_jet(n)
    result[0] = value
    return result

proc t_horner* (jet: openArray[float], h: float): float =
    for i in countdown(jet.high, jet.low):
        result = result * h + jet[i]

proc t_abs* (u: openArray[float], k: int): float =
    return if u[0] < 0.0: - u[k] else: u[k]

proc cauchy (a, b: openArray[float], k, lower, upper: int): float =
    for j in lower..upper:
        result += a[j] * b[k - j]

proc t_prod* (u, v: openArray[float], k: int): float =
    return cauchy(u, v, k, 0, k)

proc t_quot* (q: var openArray[float], u, v: openArray[float], k: int): float =
    assert v[0] != 0.0
    q[k] = (u[k] - cauchy(q, v, k, 0, k - 1)) / v[0]
    return q[k]

proc t_sqr* (u: openArray[float], k: int): float =
    result = 2.0 * cauchy(u, u, k, 0, (k - (if k mod 2 == 0: 2 else: 1)) div 2)
    if k mod 2 == 0: result += u[k div 2] * u[k div 2]

proc t_sqrt* (r: var openArray[float], u: openArray[float], k: int): float =
    assert u[0] > 0.0
    if (k == 0):
        r[k] = sqrt(u[0])
    else:
        r[k] = 2.0 * cauchy(r, r, k, 1, (k - (if k mod 2 == 0: 2 else: 1)) div 2)
        if k mod 2 == 0: r[k] += r[k div 2] * r[k div 2]
        r[k] = 0.5 * (u[k] - r[k]) / r[0]
    return r[k]

proc d_cauchy (h, u: openArray[float], k, lower, upper: int, factor: float): float =
    for j in lower..upper:
        result += factor * h[j] * float(k - j) * u[k - j] / float(k)

proc t_exp* (e: var openArray[float], u: openArray[float], k: int): float =
    if k == 0:
        e[k] = exp(u[0])
    else:
        e[k] = d_cauchy(e, u, k, 0, k - 1, 1.0)
    return e[k]

proc t_sin_cos* (s, c: var openArray[float], u: openArray[float], k: int, g: Geometry): auto =
    if k == 0:
        s[k] = if g == trig: sin(u[0]) else: sinh(u[0])
        c[k] = if g == trig: cos(u[0]) else: cosh(u[0])
    else:
        s[k] = d_cauchy(c, u, k, 0, k - 1, 1.0)
        c[k] = d_cauchy(s, u, k, 0, k - 1, if g == trig: -1.0 else: 1.0)
    return (s[k], c[k])

proc t_tan_sec2* (t, s: var openArray[float], u: openArray[float], k: int, g: Geometry): auto =
    if k == 0:
        t[k] = if g == trig: tan(u[0]) else: tanh(u[0])
        s[k] = if g == trig: sec(u[0]) * sec(u[0]) else: sech(u[0]) * sech(u[0])
    else:
        t[k] = d_cauchy(s, u, k, 0, k - 1, 1.0)
        s[k] = d_cauchy(t, t, k, 0, k - 1, if g == trig: 2.0 else: -2.0)
    return (t[k], s[k])

proc t_pow* (p: var openArray[float], u: openArray[float], a: float, k: int): float =
    assert u[0] > 0.0
    if k == 0:
        p[k] = pow(u[0], a)
    else:
        p[k] = (d_cauchy(p, u, k, 0, k - 1, a) - d_cauchy(u, p, k, 0, k - 1, 1.0)) / u[0]
    return p[k]

proc t_ln* (l: var openArray[float], u: openArray[float], k: int): float =
    assert u[0] > 0.0
    if k == 0:
        l[k] = ln(u[0])
    else:
        l[k] = (u[k] - d_cauchy(u, l, k, 1, k - 1, 1.0)) / u[0]
    return l[k]

