#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from math import sqrt, sin, cos, sinh, cosh, tan, tanh, atan, asin, acos, exp, log, asinh, acosh, atanh


def t_jet(n, value=0.0):
    jet = [0.0] * n
    jet[0] = value
    return jet

def t_horner(jet, n, h):
    result = jet[n]
    for i in range(n - 1, -1, -1):
        result = result * h + jet[i]
    return result

def t_abs(u, k):
    if k == 0:
        return abs(u[0])
    elif k == 1:
        return u[1] if u[0] > 0.0 else (- u[1] if u[0] < 0.0 else 0.0)
    return 0.0

def t_prod(u, v, k):
    return sum(u[j] * v[k - j] for j in range(k + 1))

def t_quot(q, u, v, k):
    return (u[k] - sum(v[j] * q[k - j] for j in range(1, k + 1))) / v[0]

def t_sqrt(r, u, k):
    if k == 0:
        return sqrt(u[0])
    return (u[k] - sum(r[j] * r[k - j] for j in range(1, k))) / (2.0 * r[0])

def t_pwr(p, u, a, k):
    if abs(u[0]) == 0.0:
        return 0.0
    if k == 0:
        return u[0]**a
    return sum((a * (k - j) - j) * p[j] * u[k - j] for j in range(k)) / (k * u[0])

def _ddot(v, u, k):
    return sum(j * u[j] * v[k - j] for j in range(1, k)) / k

def t_exp(e, u, k):
    if k == 0:
        return exp(u[0])
    return e[0] * u[k] + _ddot(e, u, k)

def t_ln(l, u, k):
    if k == 0:
        return log(u[0])
    return (u[k] - _ddot(u, l, k)) / u[0]

def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    sk = c[0] * u[k] + _ddot(c, u, k)
    ck = s[0] * u[k] + _ddot(s, u, k)
    return (sk, ck) if hyp else (sk, - ck)

def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyp else (tan(u[0]), 1.0 + tan(u[0])**2)
    tk = s2[0] * u[k] + _ddot(s2, u, k)
    sk = 2.0 * (t[0] * tk + _ddot(t, t, k))
    return (tk, - sk) if hyp else (tk, sk)

def t_asin(h, v, u, k, hyp=False):
    if k == 0:
        return (asinh(u[0]), cosh(asinh(u[0]))) if hyp else (asin(u[0]), cos(asin(u[0])))
    sk = (u[k] - _ddot(v, h, k)) / v[0]
    ck = u[0] * sk + _ddot(u, h, k)
    return (sk, ck) if hyp else (sk, - ck)

def t_acos(h, v, u, k, hyp=False):
    if k == 0:
        return (acosh(u[0]), sinh(acosh(u[0]))) if hyp else (acos(u[0]), - sin(acos(u[0])))
    sk = (u[k] + _ddot(v, h, k)) / v[0]
    ck = - u[0] * sk - _ddot(u, h, k)
    return (sk, ck) if hyp else (sk, - ck)

def t_atan(h, v, u, k, hyp=False):
    if k == 0:
        return (atanh(u[0]), 1.0 - u[0]**2) if hyp else (atan(u[0]), 1.0 + u[0]**2)
    tk, sk = (u[k] - _ddot(v, h, k)) / v[0], 2.0 * (u[0] * u[k] + _ddot(u, u, k))
    return (tk, - sk) if hyp else (tk, sk)


class Series:

    def __init__(self, jet):
        self.n = len(jet)
        self.jet = [jet[k] for k in range(self.n)]

    @classmethod
    def get(cls, order, value=0.0):
        return cls(t_jet(order, value))

    def __str__(self):
        return ''.join(f"{self.jet[i]:+.6e} " for i in range(self.n))

    def __abs__(self):
        return Series([t_abs(self.jet, k) for k in range(self.n)])

    def __pos__(self):
        return Series([self.jet[k] for k in range(self.n)])

    def __neg__(self):
        return Series([- self.jet[k] for k in range(self.n)])

    def __invert__(self):  # override - returns a derivative Series
        d = t_jet(self.n, self.jet[0])
        fac = 1
        for i in range(1, self.n):
            fac *= i
            d[i] = fac * self.jet[i]
        return Series(d)

    def __add__(self, other):
        if isinstance(other, Series):
            return Series([self.jet[k] + other.jet[k] for k in range(self.n)])
        add_jet = [self.jet[k] for k in range(self.n)]
        add_jet[0] += other
        return Series(add_jet)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Series):
            return Series([self.jet[k] - other.jet[k] for k in range(self.n)])
        sub_jet = [self.jet[k] for k in range(self.n)]
        sub_jet[0] -= other
        return Series(sub_jet)

    def __rsub__(self, other):
        sub_jet = [- self.jet[k] for k in range(self.n)]
        sub_jet[0] += other
        return Series(sub_jet)

    def __mul__(self, other):
        if isinstance(other, Series):
            return Series([t_prod(self.jet, other.jet, k) for k in range(self.n)])
        return Series([self.jet[k] * other for k in range(self.n)])

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Series):
            assert abs(other.val) != 0.0, f"other.val = {other.val}"
            div_jet = t_jet(self.n)
            for k in range(self.n):
                div_jet[k] = t_quot(div_jet, self.jet, other.jet, k)
            return Series(div_jet)
        assert abs(other) != 0.0, f"other = {other}"
        return Series([self.jet[k] / other for k in range(self.n)])

    def __rtruediv__(self, other):
        assert abs(self.val) != 0.0, f"self.val = {self.val}"
        other_jet = t_jet(self.n, other)
        rdiv_jet = t_jet(self.n)
        for k in range(self.n):
            rdiv_jet[k] = t_quot(rdiv_jet, other_jet, self.jet, k)
        return Series(rdiv_jet)

    def __pow__(self, other):
        assert self.val > 0.0, f"self.val = {self.val}"
        if isinstance(other, Series):
            return (self.ln * other).exp
        pow_jet = t_jet(self.n)
        for k in range(self.n):
            pow_jet[k] = t_pwr(pow_jet, self.jet, other, k)
        return Series(pow_jet)

    def __rpow__(self, other):
        assert other > 0.0, f"other = {other}"
        return (log(other) * self).exp

    def _single(self, fun):
        jet = t_jet(self.n)
        for k in range(self.n):
            jet[k] = fun(jet, self.jet, k)
        return Series(jet)

    def _double(self, fun, hyp=False):
        jet_a, jet_b = t_jet(self.n), t_jet(self.n)
        for k in range(self.n):
            jet_a[k], jet_b[k] = fun(jet_a, jet_b, self.jet, k, hyp)
        return Series(jet_a), Series(jet_b)

    @property
    def sqrt(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return self._single(t_sqrt)

    @property
    def exp(self):
        return self._single(t_exp)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return self._single(t_ln)

    @property
    def sin(self):
        return self._double(t_sin_cos)[0]

    @property
    def cos(self):
        return self._double(t_sin_cos)[1]

    @property
    def sin_cos(self):
        return self._double(t_sin_cos)

    @property
    def tan(self):
        return self._double(t_tan_sec2)[0]

    @property
    def sinh(self):
        return self._double(t_sin_cos, hyp=True)[0]

    @property
    def cosh(self):
        return self._double(t_sin_cos, hyp=True)[1]

    @property
    def sinh_cosh(self):
        return self._double(t_sin_cos, hyp=True)

    @property
    def tanh(self):
        return self._double(t_tan_sec2, hyp=True)[0]

    @property
    def asin(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._double(t_asin)[0]

    @property
    def acos(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._double(t_acos)[0]

    @property
    def atan(self):
        return self._double(t_atan)[0]

    @property
    def asinh(self):
        return self._double(t_asin, hyp=True)[0]

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return self._double(t_acos, hyp=True)[0]

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return self._double(t_atan, hyp=True)[0]

    @property
    def val(self):
        return self.jet[0]

    @property
    def var(self):
        jet = t_jet(self.n, self.val)
        jet[1] = 1.0
        return Series(jet)


class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value=0.0):
        return cls(value, 0.0)

    def __str__(self):
        return f"{self.val:+.6e} {self.der:+.6e}"

    def __abs__(self):
        return Dual(abs(self.val), self.der if self.val > 0.0 else (- self.der if self.val < 0.0 else 0.0))

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val + other.val, self.der + other.der)
        return Dual(self.val + other, self.der)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val - other.val, self.der - other.der)
        return Dual(self.val - other, self.der)

    def __rsub__(self, other):
        return Dual(- self.val + other, - self.der)

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val * other.val, self.der * other.val + self.val * other.der)
        return Dual(self.val * other, self.der * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, Dual):
            assert abs(other.val) != 0.0, f"other.val = {other.val}"
            return Dual(self.val / other.val, (self.der * other.val - self.val * other.der) / other.val**2)
        assert abs(other) != 0.0, f"other = {other}"
        return Dual(self.val / other, self.der / other)

    def __rtruediv__(self, other):
        assert abs(self.val) != 0.0, f"self.val = {self.val}"
        return Dual(other / self.val, - other * self.der / self.val**2)

    def __pow__(self, other):
        assert self.val > 0.0, f"self.val = {self.val}"
        if isinstance(other, Dual):
            return (self.ln * other).exp
        return Dual(self.val**other, other * self.val**(other - 1.0) * self.der)

    def __rpow__(self, other):
        assert other > 0.0, f"other = {other}"
        return (log(other) * self).exp

    @property
    def sqrt(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        sqrt_val = sqrt(self.val)
        return Dual(sqrt_val, self.der / (2.0 * sqrt_val))

    @property
    def exp(self):
        exp_val = exp(self.val)
        return Dual(exp_val, self.der * exp_val)

    @property
    def ln(self):
        assert self.val > 0.0, f"self.val = {self.val}"
        return Dual(log(self.val), self.der / self.val)

    @property
    def sin(self):
        return Dual(sin(self.val), self.der * cos(self.val))

    @property
    def cos(self):
        return Dual(cos(self.val), - self.der * sin(self.val))

    @property
    def tan(self):
        return Dual(tan(self.val), self.der * (1.0 + tan(self.val)**2))

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.der * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.der * sinh(self.val))

    @property
    def tanh(self):
        return Dual(tanh(self.val), self.der * (1.0 - tanh(self.val)**2))

    @property
    def asin(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(asin(self.val), self.der / sqrt(1.0 - self.val**2))

    @property
    def acos(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(acos(self.val), - self.der / sqrt(1.0 - self.val**2))

    @property
    def atan(self):
        return Dual(atan(self.val), self.der / (1.0 + self.val**2))

    @property
    def asinh(self):
        return Dual(asinh(self.val), self.der / sqrt(self.val**2 + 1.0))

    @property
    def acosh(self):
        assert self.val > 1.0, f"self.val = {self.val}"
        return Dual(acosh(self.val), self.der / sqrt(self.val**2 - 1.0))

    @property
    def atanh(self):
        assert abs(self.val) < 1.0, f"self.val = {self.val}"
        return Dual(atanh(self.val), self.der / (1.0 - self.val**2))

    @property
    def var(self):
        return Dual(self.val, 1.0)


print(__name__ + " module loaded", file=stderr)
