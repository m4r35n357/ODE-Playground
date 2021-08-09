#
#  (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr, argv
from math import fsum, sin, cos, sinh, cosh, tan, tanh, exp, log, asinh, asin, acosh, acos, atanh, atan, sqrt
from collections import namedtuple


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

def _i_cauchy(g, u, f, k, sign=True):
    return ((u[k] if sign else - u[k]) - _d_cauchy(g, f, k, 1, k - 1)) / g[0]

def t_asin(h, g, u, k, hyp=False):
    if k == 0:
        return (asinh(u[0]), sqrt(u[0]**2 + 1.0)) if hyp else (asin(u[0]), sqrt(1.0 - u[0]**2))
    h[k] = _i_cauchy(g, u, h, k)
    g[k] = _d_cauchy(u, h, k, 0, k - 1)
    return (h[k], g[k]) if hyp else (h[k], - g[k])

def t_acos(h, g, u, k, hyp=False):
    if k == 0:
        return (acosh(u[0]), sqrt(u[0]**2 - 1.0)) if hyp else (acos(u[0]), sqrt(1.0 - u[0]**2))
    h[k] = _i_cauchy(g, u, h, k, sign=hyp)
    g[k] = _d_cauchy(u, h, k, 0, k - 1)
    return h[k], g[k]

def t_atan(h, g, u, k, hyp=False):
    if k == 0:
        return (atanh(u[0]), 1.0 - u[0]**2) if hyp else (atan(u[0]), 1.0 + u[0]**2)
    h[k] = _i_cauchy(g, u, h, k)
    g[k] = 2.0 * _d_cauchy(u, u, k, 0, k - 1)
    return (h[k], - g[k]) if hyp else (h[k], g[k])


class Components(namedtuple('ParametersType', ['x', 'y', 'z'])):
    pass

def output(x, y, z, t):
    print(f'{x:+.{Context.places}e} {y:+.{Context.places}e} {z:+.{Context.places}e} {t:.5e}')

def tsm(ode, get_p, get_i):
    Context.places, n, δt, n_steps = int(argv[1]), int(argv[2]), float(argv[3]), int(argv[4])  # controls
    x, y, z = t_jet(n + 1), t_jet(n + 1), t_jet(n + 1)  # coordinate jets
    x[0], y[0], z[0] = float(argv[6]), float(argv[7]), float(argv[8])  # initial values
    steps = range(1, n_steps + 1)
    index = range(n)
    p = get_p(n)
    i = None if get_i is None else get_i(n)
    output(x[0], y[0], z[0], 0.0)
    for step in steps:
        for k in index:
            c = ode(x, y, z, p, i, k)
            x[k + 1] = c.x / (k + 1)
            y[k + 1] = c.y / (k + 1)
            z[k + 1] = c.z / (k + 1)
        x[0], y[0], z[0] = t_horner(x, δt), t_horner(y, δt), t_horner(z, δt)
        output(x[0], y[0], z[0], step * δt)


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
            jet[k] = t_inv(jet, self.jet, k) * o
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

    def _trig_hyp(self, f, hyp=False):
        a, b = t_jet(self.n), t_jet(self.n)
        for k in self.index:
            a[k], b[k] = f(a, b, self.jet, k, hyp)
        return Series(a), Series(b)

    @property
    def exp(self):
        return self._exp_log_sqrt(t_exp)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return self._exp_log_sqrt(t_ln)

    @property
    def sin(self):
        return self._trig_hyp(t_sin_cos)[0]

    @property
    def cos(self):
        return self._trig_hyp(t_sin_cos)[1]

    @property
    def tan(self):
        return self._trig_hyp(t_tan_sec2)[0]

    @property
    def sec2(self):
        return self._trig_hyp(t_tan_sec2)[1]

    @property
    def sinh(self):
        return self._trig_hyp(t_sin_cos, hyp=True)[0]

    @property
    def cosh(self):
        return self._trig_hyp(t_sin_cos, hyp=True)[1]

    @property
    def tanh(self):
        return self._trig_hyp(t_tan_sec2, hyp=True)[0]

    @property
    def sech2(self):
        return self._trig_hyp(t_tan_sec2, hyp=True)[1]

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
        return self._trig_hyp(t_asin, hyp=True)[0]

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_acos, hyp=True)[0]

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._trig_hyp(t_atan, hyp=True)[0]

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

    def __init__(self, value, derivative):
        self.val = value
        self.dot = derivative

    @classmethod
    def get(cls, value=0.0):
        return cls(value if isinstance(value, float) else float(value), 0.0)

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
            return Dual(self.val ** o, self.dot * o * self.val ** (o - 1))
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
