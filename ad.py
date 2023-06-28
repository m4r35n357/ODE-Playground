#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from math import fsum, sin, cos, sinh, cosh, tan, tanh, exp, log, asinh, asin, acosh, acos, atanh, atan, sqrt
from collections import namedtuple
from time import clock_gettime, CLOCK_MONOTONIC


def t_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

def t_horner(jet, h):
    result = 0.0
    for term in reversed(jet):
        result = result * h + term
    return result

def t_const(a, k):
    return a if k == 0 else 0.0

def t_abs(u, k):
    return - u[k] if u[0] < 0.0 else u[k]

def fa(a, b, k0, k1, k):
    return fsum(a[j] * b[k - j] for j in range(k0, k1))

def t_prod(u, v, k):
    return fa(u, v, 0, k + 1, k)

def t_quot(q, u, v, k):
    return ((u[0] if u else 1.0) if k == 0 else (u[k] if u else 0.0) - fa(q, v, 0, k, k)) / v[0]

def half(k):
    return 1 + (k - (1 if k % 2 else 2)) // 2

def rem(a, k):
    return 0.0 if k % 2 else a[k // 2] * a[k // 2]

def t_sqr(u, k):
    return 2.0 * fa(u, u, 0, half(k), k) + rem(u, k)

def t_sqrt(r, u, k):
    return sqrt(u[k]) if k == 0 else 0.5 * (u[k] - 2.0 * fa(r, r, 1, half(k), k) - rem(r, k)) / r[0]

def fb(df_du, u, k):
    return fsum(df_du[j] * (k - j) * u[k - j] for j in range(k)) / k

def t_exp(e, u, k):
    return exp(u[k]) if k == 0 else fb(e, u, k)

def t_sin_cos(s, c, u, k, trig=True):
    if k == 0:
        return (sin(u[k]), cos(u[k])) if trig else (sinh(u[k]), cosh(u[k]))
    s[k] = fb(c, u, k)
    c[k] = fb(s, u, k)
    return s[k], -c[k] if trig else c[k]

def t_tan_sec2(t, s2, u, k, trig=True):
    if k == 0:
        t[k] = tan(u[k]) if trig else tanh(u[k])
        return (t[k], 1.0 + t[k] * t[k]) if trig else (t[k], 1.0 - t[k] * t[k])
    t[k] = fb(s2, u, k)
    s2[k] = fb(t, t, k)
    return t[k], (2.0 if trig else -2.0) * s2[k]

def t_pwr(p, u, a, k):
    return u[k]**a if k == 0 else fsum((a * (k - j) - j) * p[j] * u[k - j] for j in range(k)) / (k * u[0])

def fc(f, du_df, u, k, flag=False):
    return (u[k] + (1.0 if flag else -1.0) * fsum(du_df[j] * (k - j) * f[k - j] for j in range(1, k)) / k) / du_df[0]

def t_ln(ln, u, k):
    return log(u[k]) if k == 0 else fc(ln, u, u, k)

def t_asin(a, g, u, k, trig=True):
    if k == 0:
        return (asin(u[k]), sqrt(1.0 - u[k] * u[k])) if trig else (asinh(u[k]), sqrt(1.0 + u[k] * u[k]))
    a[k] = fc(a, g, u, k)
    g[k] = fb(u, a, k)
    return a[k], -g[k] if trig else g[k]

def t_acos(a, g, u, k, trig=True):
    if k == 0:
        return (acos(u[k]), - sqrt(1.0 - u[k] * u[k])) if trig else (acosh(u[k]), sqrt(u[k] * u[k] - 1.0))
    a[k] = fc(a, g, u, k, trig)
    g[k] = fb(u, a, k)
    return a[k], g[k]

def t_atan(a, g, u, k, trig=True):
    if k == 0:
        return (atan(u[k]), 1.0 + u[k] * u[k]) if trig else (atanh(u[k]), 1.0 - u[k] * u[k])
    a[k] = fc(a, g, u, k)
    g[k] = 2.0 * fb(u, u, k)
    return a[k], g[k] if trig else -g[k]


def t_out(dp, x, y, z, t, cpu):
    print(f'{x:+.{dp}e} {y:+.{dp}e} {z:+.{dp}e} {t:.5e} _ _ _ {cpu:.5e}')


def tsm(ode, places, n, h, steps, x0, y0, z0, p):
    x = t_jet(n + 1, x0)
    y = t_jet(n + 1, y0)
    z = t_jet(n + 1, z0)
    t0 = clock_gettime(CLOCK_MONOTONIC)
    for step in range(steps):
        for k in range(n):
            c = ode(x, y, z, p, k)
            x[k + 1] = c.x / (k + 1)
            y[k + 1] = c.y / (k + 1)
            z[k + 1] = c.z / (k + 1)
        t_out(places, x[0], y[0], z[0], step * h, clock_gettime(CLOCK_MONOTONIC) - t0)
        x[0] = t_horner(x, h)
        y[0] = t_horner(y, h)
        z[0] = t_horner(z, h)
    t_out(places, x[0], y[0], z[0], steps * h, clock_gettime(CLOCK_MONOTONIC) - t0)


def rk4(ode, places, skip, h, steps, x, y, z, p):
    t0 = clock_gettime(CLOCK_MONOTONIC)
    for step in range(steps):
        k1 = ode(x, y, z, p)
        k2 = ode(x + 0.5 * k1.x * h, y + 0.5 * k1.y * h, z + 0.5 * k1.z * h, p)
        k3 = ode(x + 0.5 * k2.x * h, y + 0.5 * k2.y * h, z + 0.5 * k2.z * h, p)
        k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p)
        x += h * (k1.x + 2.0 * (k2.x + k3.x) + k4.x) / 6.0
        y += h * (k1.y + 2.0 * (k2.y + k3.y) + k4.y) / 6.0
        z += h * (k1.z + 2.0 * (k2.z + k3.z) + k4.z) / 6.0
        if step % skip == 0:
            t_out(places, x, y, z, step * h, clock_gettime(CLOCK_MONOTONIC) - t0)
    t_out(places, x, y, z, steps * h, clock_gettime(CLOCK_MONOTONIC) - t0)


class Components(namedtuple('ParametersType', ['x', 'y', 'z'])):
    pass


class Context:
    places = 3


class Series:

    def __init__(self, jet):
        self.n = len(jet)
        self.index = range(self.n)
        self.jet = jet[:]

    @classmethod
    def get(cls, size, value=0.0):
        return cls(t_jet(size, value))

    def __str__(self):
        return ''.join(f'{term:+.{Context.places}e} ' for term in self.jet)

    def __abs__(self):
        return - self if self.val < 0.0 else + self

    def __pos__(self):
        return Series(self.jet[:])

    def __neg__(self):
        return Series([- term for term in self.jet])

    def __invert__(self):  # override - returns a derivative Series
        derivatives = t_jet(self.n, self.val)
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
        return self + (- o)

    def __rsub__(self, o):
        return - self + o

    def __mul__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return Series([t_prod(self.jet, o.jet, k) for k in self.index])
        elif isinstance(o, (float, int)):
            return Series([term * o for term in self.jet])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self * o

    def __truediv__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            assert abs(o.val) != 0.0, f"other.val = {o.val}"
            jet = t_jet(self.n)
            for k in self.index:
                jet[k] = t_quot(jet, self.jet, o.jet, k)
            return Series(jet)
        elif isinstance(o, (float, int)):
            assert abs(o) != 0.0, f"other = {o}"
            return Series([term / o for term in self.jet])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        assert self.val != 0.0, f"self.val = {self.val}"
        jet = t_jet(self.n)
        for k in self.index:
            # noinspection PyTypeChecker
            jet[k] = t_quot(jet, None, self.jet, k) * o
        return Series(jet)

    def __pow__(self, o):
        if isinstance(o, int):
            i_pow = Series(self.jet)
            for _ in range(abs(o) - 1):
                i_pow = i_pow * self
            return i_pow if o > 0 else (1.0 / i_pow if o < 0 else Series(t_jet(self.n, 1.0)))
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
            jet = t_jet(self.n)
            for k in self.index:
                jet[k] = t_pwr(jet, self.jet, o, k)
            return Series(jet)
        elif isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return (self.ln * o).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self * log(o)).exp

    def _exp_log_sqrt(self, f):
        a = t_jet(self.n)
        for k in self.index:
            a[k] = f(a, self.jet, k)
        return Series(a)

    def _trig_hyp(self, f, trig=True):
        a, b = t_jet(self.n), t_jet(self.n)
        for k in self.index:
            a[k], b[k] = f(a, b, self.jet, k, trig)
        return Series(a), Series(b)

    @property
    def exp(self):
        return self._exp_log_sqrt(t_exp)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return self._exp_log_sqrt(t_ln)

    @property
    def sin_cos(self):
        return self._trig_hyp(t_sin_cos)

    @property
    def sin(self):
        return self._trig_hyp(t_sin_cos)[0]

    @property
    def cos(self):
        return self._trig_hyp(t_sin_cos)[1]

    @property
    def tan_sec2(self):
        return self._trig_hyp(t_tan_sec2)

    @property
    def tan(self):
        return self._trig_hyp(t_tan_sec2)[0]

    @property
    def sec2(self):
        return self._trig_hyp(t_tan_sec2)[1]

    @property
    def sinh_cosh(self):
        return self._trig_hyp(t_sin_cos, trig=False)

    @property
    def sinh(self):
        return self._trig_hyp(t_sin_cos, trig=False)[0]

    @property
    def cosh(self):
        return self._trig_hyp(t_sin_cos, trig=False)[1]

    @property
    def tanh_sech2(self):
        return self._trig_hyp(t_tan_sec2, trig=False)

    @property
    def tanh(self):
        return self._trig_hyp(t_tan_sec2, trig=False)[0]

    @property
    def sech2(self):
        return self._trig_hyp(t_tan_sec2, trig=False)[1]

    @property
    def asin(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_asin)[0]

    @property
    def acos(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_acos)[0]

    @property
    def atan(self):
        return self._trig_hyp(t_atan)[0]

    @property
    def asinh(self):
        return self._trig_hyp(t_asin, trig=False)[0]

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_acos, trig=False)[0]

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_atan, trig=False)[0]

    @property
    def sqr(self):
        return Series([t_sqr(self.jet, k) for k in self.index])

    @property
    def sqrt(self):
        assert abs(self.val) > 0.0, f"self.val = {self.val}"
        return self._exp_log_sqrt(t_sqrt)

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
        return f'{self.val:+.{Context.places}e} {self.dot:+.{Context.places}e}'

    def __abs__(self):
        return - self if self.val < 0.0 else + self

    def __pos__(self):
        return Dual(self.val, self.dot)

    def __neg__(self):
        return Dual(- self.val, - self.dot)

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val + o.val, self.dot + o.dot)
        elif isinstance(o, (float, int)):
            return Dual(self.val + o, self.dot)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __radd__(self, o):
        return self + o

    def __sub__(self, o):
        return self + (- o)

    def __rsub__(self, o):
        return - self + o

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
            assert abs(o.val) != 0.0, f"other.val = {o.val}"
            return Dual(self.val / o.val, (self.dot * o.val - self.val * o.dot) / o.val ** 2)
        elif isinstance(o, (float, int)):
            assert abs(o) != 0.0, f"other = {o}"
            return Dual(self.val / o, self.dot / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        assert self.val != 0.0, f"self.val = {self.val}"
        return Dual(o / self.val, - self.dot * o / self.val ** 2)

    def __pow__(self, o):
        if isinstance(o, int):
            i_pow = Dual(self.val, self.dot)
            for _ in range(abs(o) - 1):
                i_pow = i_pow * self
            return i_pow if o > 0 else (1.0 / i_pow if o < 0 else Dual(1.0, 0.0))
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
            pwr = self.val ** o
            return Dual(pwr, self.dot * o * pwr / self.val)
        elif isinstance(o, Dual):
            return (self.ln * o).exp
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self * log(o)).exp

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
        return Dual(t, self.dot * (1.0 + t ** 2))

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.dot * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.dot * sinh(self.val))

    @property
    def tanh(self):
        t = tanh(self.val)
        return Dual(t, self.dot * (1.0 - t ** 2))

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


print(f'{__name__} module loaded', file=stderr)
