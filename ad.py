#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from math import fsum, sin, cos, sinh, cosh, tan, tanh, exp, log, asinh, asin, acosh, acos, atanh, atan, sqrt
from sys import stderr


def t_jet(n, value=0.0):
    return [value if isinstance(value, float) else float(value)] + [0.0] * (n - 1)

def t_horner(jet, h):
    result = 0.0
    for term in reversed(jet):
        result = result * h + term
    return result

def t_abs(u, k):
    return (1.0 if u[0] > 0.0 else (- 1.0 if u[0] < 0.0 else 0.0)) * u[k]

def t_prod(u, v, k):
    return fsum(u[j] * v[k - j] for j in range(k + 1))

def t_quot(q, u, v, k):
    return (u[k] - fsum(q[j] * v[k - j] for j in range(k))) / v[0]

def t_sqr(u, k):
    f_sum = 2.0 * fsum(u[j] * u[k - j] for j in (range((k - 1) // 2 + 1) if k % 2 == 1 else range((k - 2) // 2 + 1)))
    return f_sum if k % 2 == 1 else f_sum + u[k // 2]**2

def t_sqrt(r, u, k):
    if k == 0:
        return sqrt(u[0])
    f_sum = 2.0 * fsum(r[j] * r[k - j] for j in (range(1, (k - 1) // 2 + 1) if k % 2 == 1 else range(1, (k - 2) // 2 + 1)))
    return 0.5 * (u[k] - f_sum if k % 2 == 1 else u[k] - f_sum - r[k // 2]**2) / r[0]

def t_pwr(p, u, a, k):
    return u[0]**a if k == 0 else fsum((a * (k - j) - j) * p[j] * u[k - j] for j in range(k)) / (k * u[0])

def t_exp(e, u, k):
    return exp(u[0]) if k == 0 else fsum((k - j) * e[j] * u[k - j] for j in range(k)) / k

def t_ln(l, u, k):
    return log(u[0]) if k == 0 else (u[k] - fsum(u[j] * (k - j) * l[k - j] for j in range(1, k)) / k) / u[0]

def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    s[k] = fsum((k - j) * c[j] * u[k - j] for j in range(k)) / k
    c[k] = fsum((k - j) * s[j] * u[k - j] for j in range(k)) / k
    return (s[k], c[k]) if hyp else (s[k], - c[k])

def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyp else (tan(u[0]), 1.0 + tan(u[0])**2)
    t[k] = fsum((k - j) * s2[j] * u[k - j] for j in range(k)) / k
    s2[k] = 2.0 * fsum((k - j) * t[j] * t[k - j] for j in range(k)) / k
    return (t[k], - s2[k]) if hyp else (t[k], s2[k])

def t_asin(h, v, u, k, hyp=False):
    if k == 0:
        return (asinh(u[0]), sqrt(u[0]**2 + 1.0)) if hyp else (asin(u[0]), sqrt(1.0 - u[0]**2))
    h[k] = (u[k] - fsum(j * h[j] * v[k - j] for j in range(1, k)) / k) / v[0]
    v[k] = u[0] * h[k] + fsum(j * h[j] * u[k - j] for j in range(1, k)) / k
    return (h[k], v[k]) if hyp else (h[k], - v[k])

def t_acos(h, v, u, k, hyp=False):
    if k == 0:
        return (acosh(u[0]), sqrt(u[0]**2 - 1.0)) if hyp else (acos(u[0]), - sqrt(1.0 - u[0]**2))
    h[k] = (u[k] + fsum(j * h[j] * v[k - j] for j in range(1, k)) / k) / v[0]
    v[k] = u[0] * h[k] + fsum(j * h[j] * u[k - j] for j in range(1, k)) / k
    return (h[k], - v[k]) if hyp else (h[k], v[k])

def t_atan(h, v, u, k, hyp=False):
    if k == 0:
        return (atanh(u[0]), 1.0 - u[0]**2) if hyp else (atan(u[0]), 1.0 + u[0]**2)
    h[k] = (u[k] - fsum(j * h[j] * v[k - j] for j in range(1, k)) / k) / v[0]
    v[k] = 2.0 * (u[0] * u[k] + fsum(j * u[j] * u[k - j] for j in range(1, k)) / k)
    return (h[k], - v[k]) if hyp else (h[k], v[k])


class Series:

    def __init__(self, jet):
        self.n = len(jet)
        self.index = range(self.n)
        self.jet = jet[:]

    @classmethod
    def get(cls, order, value=0.0):
        return cls(t_jet(order, value))

    def __str__(self):
        return ''.join(f'{term:+.{Context.places}e} ' for term in self.jet)

    def __abs__(self):
        return Series([t_abs(self.jet, k) for k in self.index])

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
        return self.__add__(o)

    def __sub__(self, o):
        return self.__add__(- o)

    def __rsub__(self, o):
        return (- self).__add__(o)

    def __mul__(self, o):
        if isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return Series([t_prod(self.jet, o.jet, k) for k in self.index])
        elif isinstance(o, (float, int)):
            return Series([term * o for term in self.jet])
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self.__mul__(o)

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
            return self.__mul__(1.0 / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        return Series(t_jet(self.n, o)).__truediv__(self)

    def __pow__(self, o):
        if isinstance(o, int):
            if o == 0:
                return Series(t_jet(self.n, 1.0))
            i_pow = Series(self.jet) if o > 0 else Series(self.jet).__rtruediv__(1.0)  # pragma: no mutate
            multiplier = Series(i_pow.jet)
            for _ in range(abs(o) - 1):
                i_pow = i_pow.__mul__(multiplier)
            return i_pow
        elif isinstance(o, Series):
            assert o.n == self.n, f"Size mismatch - self: {self.n}, other: {o.n}"
            return (self.ln.__mul__(o)).exp
        elif isinstance(o, float):
            assert self.val > 0.0, f"self.val = {self.val}"
            jet = t_jet(self.n)
            for k in self.index:
                jet[k] = t_pwr(jet, self.jet, o, k)
            return Series(jet)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self.__mul__(log(o))).exp

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
        assert self.n > 1, f"Order-1 series {self.val} cannot be a variable"
        return Series([self.val] + [1.0] + [0.0] * (self.n - 2))


class Context:
    places = 3


print(f'{__name__} module loaded', file=stderr)
