#
# Taylor Series Method for solving ODEs, with a Taylor Arithmetic class & a Dual Numbers class
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from math import sin, cos, sinh, cosh, tan, tanh, exp, log, asinh, asin, acosh, acos, atanh, atan, sqrt, fsum
from collections import namedtuple
from time import clock_gettime, CLOCK_MONOTONIC

def tsm_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

def horner(jet, h):
    result = 0.0
    for term in reversed(jet):
        result = result * h + term
    return result

def _out_(dp, x, y, z, t, cpu):
    print(f'{x:+.{dp}e} {y:+.{dp}e} {z:+.{dp}e} {t:.5e} _ _ _ {cpu:.5e}')

def tsm(ode, places, n, h, steps, x0, y0, z0, p):
    x = tsm_jet(n + 1, x0)
    y = tsm_jet(n + 1, y0)
    z = tsm_jet(n + 1, z0)
    t0 = clock_gettime(CLOCK_MONOTONIC)
    for step in range(steps):
        _out_(places, x[0], y[0], z[0], step * h, clock_gettime(CLOCK_MONOTONIC) - t0)
        for k in range(n):
            v = ode(x, y, z, p, k)
            x[k + 1] = v.x / (k + 1)
            y[k + 1] = v.y / (k + 1)
            z[k + 1] = v.z / (k + 1)
        x[0] = horner(x, h)
        y[0] = horner(y, h)
        z[0] = horner(z, h)
    _out_(places, x[0], y[0], z[0], steps * h, clock_gettime(CLOCK_MONOTONIC) - t0)

def rk4(ode, places, h, steps, x, y, z, p):
    t0 = clock_gettime(CLOCK_MONOTONIC)
    for step in range(steps):
        _out_(places, x, y, z, step * h, clock_gettime(CLOCK_MONOTONIC) - t0)
        k1 = ode(x, y, z, p)
        k2 = ode(x + 0.5 * k1.x * h, y + 0.5 * k1.y * h, z + 0.5 * k1.z * h, p)
        k3 = ode(x + 0.5 * k2.x * h, y + 0.5 * k2.y * h, z + 0.5 * k2.z * h, p)
        k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p)
        x += h * (k1.x + 2.0 * (k2.x + k3.x) + k4.x) / 6.0
        y += h * (k1.y + 2.0 * (k2.y + k3.y) + k4.y) / 6.0
        z += h * (k1.z + 2.0 * (k2.z + k3.z) + k4.z) / 6.0
    _out_(places, x, y, z, steps * h, clock_gettime(CLOCK_MONOTONIC) - t0)

def t_const(v, k):
    return v if k == 0 else 0.0

def t_abs(u, k):
    return - u[k] if u[0] < 0.0 else u[k]

def _cauchy_(b, a, k, k0, k1):
    return fsum(b[j] * a[k - j] for j in range(k0, k1 + 1))

def _half_(a, k, k0):
    return 2.0 * _cauchy_(a, a, k, k0, (k - (1 if k % 2 else 2)) // 2) + (0.0 if k % 2 else a[k // 2] * a[k // 2])

def _chain_(b, a, k, k0):
    return fsum(b[j] * (k - j) * a[k - j] for j in range(k0, k)) / k

def _fk_(df_du, u, k, scale=1.0):
    return _chain_(df_du, u, k, 0) * scale

def _uk_(df_du, u, k, f_k, scale=1.0):
    return (f_k - _chain_(df_du, u, k, 1) * scale) / df_du[0]

def t_mul(u, v, k):
    return _cauchy_(u, v, k, 0, k)

def t_div(q, u, v, k):
    assert v[0] != 0.0
    q[k] = u[k] if k == 0 else u[k] - _cauchy_(q, v, k, 0, k - 1)
    return q[k] / v[0]

def t_sqr(u, k):
    return _half_(u, k, 0)

def t_sqrt(r, u, k):
    assert u[0] > 0.0
    r[k] = sqrt(u[k]) if k == 0 else 0.5 * (u[k] - _half_(r, k, 1)) / r[0]
    return r[k]

def t_exp(e, u, k):
    e[k] = exp(u[k]) if k == 0 else _fk_(e, u, k)
    return e[k]

def t_sin_cos(s, c, u, k, trig=True):
    if k == 0:
        s[k] = sin(u[k]) if trig else sinh(u[k])
        c[k] = cos(u[k]) if trig else cosh(u[k])
    else:
        s[k] = _fk_(c, u, k)
        c[k] = _fk_(s, u, k, -1.0 if trig else 1.0)
    return s[k], c[k]

def t_tan_sec2(t, s, u, k, trig=True):
    if k == 0:
        t[k] = tan(u[k]) if trig else tanh(u[k])
        s[k] = 1.0 + t[k] * t[k] if trig else 1.0 - t[k] * t[k]
    else:
        t[k] = _fk_(s, u, k)
        s[k] = _fk_(t, t, k, 2.0 if trig else -2.0)
    return t[k], s[k]

def t_ln(u, e, k):
    assert e[0] > 0.0
    u[k] = log(e[k]) if k == 0 else _uk_(e, u, k, e[k])
    return u[k]

def t_asin(u, c, s, k, trig=True):
    assert -1.0 < s[0] < 1.0 if trig else True
    if k == 0:
        u[k] = asin(s[k]) if trig else asinh(s[k])
        c[k] = cos(u[k]) if trig else cosh(u[k])
    else:
        u[k] = _uk_(c, u, k, s[k])
        c[k] = _fk_(s, u, k, -1.0 if trig else 1.0)
    return u[k], c[k]

def t_acos(u, s, c, k, trig=True):
    assert -1.0 < c[0] < 1.0 if trig else c[0] > 1.0
    if k == 0:
        u[k] = acos(c[k]) if trig else acosh(c[k])
        s[k] = -sin(u[k]) if trig else sinh(u[k])
    else:
        u[k] = _uk_(s, u, k, c[k], -1.0 if trig else 1.0)
        s[k] = _fk_(c, u, k)
    return u[k], s[k]

def t_atan(u, s, t, k, trig=True):
    assert True if trig else -1.0 < t[0] < 1.0
    if k == 0:
        u[k] = atan(t[k]) if trig else atanh(t[k])
        s[k] = 1.0 + t[k] * t[k] if trig else 1.0 - t[k] * t[k]
    else:
        u[k] = _uk_(s, u, k, t[k])
        s[k] = _fk_(t, t, k, 2.0 if trig else -2.0)
    return u[k], s[k]

def t_pwr(p, u, a, k):
    assert u[0] > 0.0
    p[k] = u[k]**a if k == 0 else _uk_(u, p, k, _fk_(p, u, k, a))
    return p[k]


class Series:

    def __init__(self, jet):
        self.n = len(jet)
        self.index = range(self.n)
        self.jet = jet[:]

    @classmethod
    def get(cls, size, value=0.0):
        return cls(tsm_jet(size, value))

    def __str__(self):
        return ''.join(f'{term: .{Context.places}e} ' for term in self.jet)

    def __abs__(self):
        return (-self) if self.val < 0.0 else (+self)

    def __pos__(self):
        return Series(self.jet[:])

    def __neg__(self):
        return Series([-term for term in self.jet])

    def __invert__(self):  # override - returns a derivative Series
        derivatives = tsm_jet(self.n, self.val)
        factorial = 1
        for i in range(1, self.n):
            factorial *= i
            derivatives[i] = factorial * self.jet[i]
        return Series(derivatives)

    def __add__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return Series([s + o for s, o in zip(self.jet, o.jet)])
        elif isinstance(o, (float, int)):
            return Series([self.val + o] + self.jet[1:])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __radd__(self, o):
        return self + o

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return o + (-self)

    def __mul__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return Series([t_mul(self.jet, o.jet, k) for k in self.index])
        elif isinstance(o, (float, int)):
            return Series([term * o for term in self.jet])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self * o

    def __truediv__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            assert o.val != 0.0, f"other.val = {o.val}"
            jet = tsm_jet(self.n)
            for k in self.index:
                jet[k] = t_div(jet, self.jet, o.jet, k)
            return Series(jet)
        elif isinstance(o, (float, int)):
            assert o != 0.0, f"other = {o}"
            return Series([term / o for term in self.jet])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            assert self.val != 0.0, f"self.val = {self.val}"
            jet = tsm_jet(self.n)
            for k in self.index:
                jet[k] = t_div(jet, o.jet, self.jet, k)
            return Series(jet)
        elif isinstance(o, (float, int)):
            assert self.val != 0.0, f"self.val = {self.val}"
            j_o = tsm_jet(self.n, o)
            jet = tsm_jet(self.n)
            for k in self.index:
                # noinspection PyTypeChecker
                jet[k] = t_div(jet, j_o, self.jet, k)
            return Series(jet)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __pow__(self, o):
        if isinstance(o, int):
            if o == 0:
                return Series(tsm_jet(self.n, 1.0))
            i_pow = Series(self.jet)
            for _ in range(abs(o) - 1):
                i_pow *= self
            return i_pow if o > 0 else 1.0 / i_pow
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
            jet = tsm_jet(self.n)
            for k in self.index:
                jet[k] = t_pwr(jet, self.jet, o, k)
            return Series(jet)
        elif isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return (self.ln * o).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        if isinstance(o, (float, int)):
            assert o > 0.0, f"other = {o}"
            return (log(o) * self).exp
        elif isinstance(o, Series):
            assert o.val > 0.0, f"other = {o}"
            return (o.ln * self).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def _single(self, f):
        a = tsm_jet(self.n)
        for k in self.index:
            a[k] = f(a, self.jet, k)
        return Series(a)

    def _pair(self, f, trig=True):
        a, b = tsm_jet(self.n), tsm_jet(self.n)
        for k in self.index:
            a[k], b[k] = f(a, b, self.jet, k, trig)
        return Series(a), Series(b)

    @property
    def exp(self):
        return self._single(t_exp)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return self._single(t_ln)

    @property
    def sin_cos(self):
        return self._pair(t_sin_cos)

    @property
    def sin(self):
        return self.sin_cos[0]

    @property
    def cos(self):
        return self.sin_cos[1]

    @property
    def tan_sec2(self):
        return self._pair(t_tan_sec2)

    @property
    def tan(self):
        return self.tan_sec2[0]

    @property
    def sec2(self):
        return self.tan_sec2[1]

    @property
    def sinh_cosh(self):
        return self._pair(t_sin_cos, trig=False)

    @property
    def sinh(self):
        return self.sinh_cosh[0]

    @property
    def cosh(self):
        return self.sinh_cosh[1]

    @property
    def tanh_sech2(self):
        return self._pair(t_tan_sec2, trig=False)

    @property
    def tanh(self):
        return self.tanh_sech2[0]

    @property
    def sech2(self):
        return self.tanh_sech2[1]

    @property
    def asin(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._pair(t_asin)[0]

    @property
    def acos(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._pair(t_acos)[0]

    @property
    def atan(self):
        return self._pair(t_atan)[0]

    @property
    def asinh(self):
        return self._pair(t_asin, trig=False)[0]

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return self._pair(t_acos, trig=False)[0]

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._pair(t_atan, trig=False)[0]

    @property
    def sqr(self):
        return Series([t_sqr(self.jet, k) for k in self.index])

    @property
    def sqrt(self):
        assert abs(self.val) > 0.0, f"self.val = {self.val}"
        return self._single(t_sqrt)

    @property
    def val(self):
        return self.jet[0]

    @property
    def var(self):
        assert self.n > 1, f"Single-element series {self.val} cannot be a variable"
        return Series([self.val] + [1.0] + [0.0] * (self.n - 2))


class Dual:

    def __init__(self, value=0.0, derivative=0.0):
        self.val = value
        self.dot = derivative

    def __str__(self):
        return f'{self.val: .{Context.places}e} {self.dot: .{Context.places}e}'

    def __abs__(self):
        return (-self) if self.val < 0.0 else (+self)

    def __pos__(self):
        return Dual(self.val, self.dot)

    def __neg__(self):
        return Dual(-self.val, -self.dot)

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val + o.val, self.dot + o.dot)
        elif isinstance(o, (float, int)):
            return Dual(self.val + o, self.dot)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __radd__(self, o):
        return self + o

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return o + (-self)

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val * o.val, self.dot * o.val + self.val * o.dot)
        elif isinstance(o, (float, int)):
            return Dual(self.val * o, self.dot * o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self * o

    def __truediv__(self, o):
        if isinstance(o, Dual):
            assert o.val != 0.0, f"other.val = {o.val}"
            return Dual(self.val / o.val, (self.dot * o.val - self.val * o.dot) / o.val**2)
        elif isinstance(o, (float, int)):
            assert o != 0.0, f"other = {o}"
            return Dual(self.val / o, self.dot / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        if isinstance(o, Dual):
            assert self.val != 0.0, f"self.val = {self.val}"
            return Dual(o.val / self.val, (o.dot * self.val - o.val * self.dot) / self.val**2)
        elif isinstance(o, (float, int)):
            assert self.val != 0.0, f"self.val = {self.val}"
            return Dual(o / self.val, - self.dot * o / self.val**2)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __pow__(self, o):
        if isinstance(o, int):
            if o == 0:
                return Dual(1.0, 0.0)
            i_pow = Dual(self.val, self.dot)
            for _ in range(abs(o) - 1):
                i_pow *= self
            return i_pow if o > 0 else 1.0 / i_pow
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
            pwr = self.val**o
            return Dual(pwr, self.dot * o * pwr / self.val)
        elif isinstance(o, Dual):
            return (self.ln * o).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        if isinstance(o, (float, int)):
            assert o > 0.0, f"other = {o}"
            return (log(o) * self).exp
        elif isinstance(o, Dual):
            assert o.val > 0.0, f"other = {o}"
            return (o.ln * self).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    @property
    def exp(self):
        exp_val = exp(self.val)
        return Dual(exp_val, self.dot * exp_val)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return Dual(log(self.val), self.dot / self.val)

    @property
    def sin(self):
        return Dual(sin(self.val), self.dot * cos(self.val))

    @property
    def cos(self):
        return Dual(cos(self.val), - self.dot * sin(self.val))

    @property
    def tan(self):
        t = tan(self.val)
        return Dual(t, self.dot * (1.0 + t**2))

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.dot * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.dot * sinh(self.val))

    @property
    def tanh(self):
        t = tanh(self.val)
        return Dual(t, self.dot * (1.0 - t**2))

    @property
    def asin(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(asin(self.val), self.dot / sqrt(1.0 - self.val * self.val))

    @property
    def acos(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(acos(self.val), - self.dot / sqrt(1.0 - self.val * self.val))

    @property
    def atan(self):
        return Dual(atan(self.val), self.dot / (1.0 + self.val * self.val))

    @property
    def asinh(self):
        return Dual(asinh(self.val), self.dot / sqrt(self.val * self.val + 1.0))

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return Dual(acosh(self.val), self.dot / sqrt(self.val * self.val - 1.0))

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(atanh(self.val), self.dot / (1.0 - self.val * self.val))

    @property
    def sqr(self):
        return self * self

    @property
    def sqrt(self):
        assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
        r = sqrt(self.val)
        return Dual(r, self.dot * 0.5 / r)

    @property
    def var(self):
        return Dual(self.val, 1.0)


class Components(namedtuple('ParametersType', ['x', 'y', 'z'])):
    pass


class Context:
    places = 3


print(f'{__name__} module loaded', file=stderr)
