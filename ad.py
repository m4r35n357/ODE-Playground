#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from math import sin, cos, sinh, cosh, tan, tanh, exp, log

def t_jet(n, value=0.0):
    jet = [0.0] * n
    jet[0] = value if isinstance(value, float) else float(value)
    return jet

def t_horner(jet, h):
    result = jet[-1]
    for i in range(len(jet) - 2, -1, -1):
        result = result * h + jet[i]
    return result

def t_abs(u, k):
    return abs(u[0]) if k == 0 else ((u[1] if u[0] > 0.0 else (- u[1] if u[0] < 0.0 else 0.0)) if k == 1 else 0.0)

def t_prod(u, v, k):
    return sum(u[j] * v[k - j] for j in range(k + 1))

def t_quot(q, u, v, k):
    return (u[k] - sum(v[j] * q[k - j] for j in range(1, k + 1))) / v[0]

def t_pwr(p, u, a, k):
    return u[0]**a if k == 0 else sum((a * (k - j) - j) * p[j] * u[k - j] for j in range(k)) / (k * u[0])

def t_exp(e, u, k):
    return exp(u[0]) if k == 0 else sum((k - j) * e[j] * u[k - j] for j in range(k)) / k

def t_ln(l, u, k):
    return log(u[0]) if k == 0 else (u[k] - sum(j * l[j] * u[k - j] for j in range(1, k)) / k) / u[0]

def t_sin_cos(s, c, u, k, hyp=False):
    if k == 0:
        return (sinh(u[0]), cosh(u[0])) if hyp else (sin(u[0]), cos(u[0]))
    sk = sum((k - j) * c[j] * u[k - j] for j in range(k)) / k
    ck = sum((k - j) * s[j] * u[k - j] for j in range(k)) / k
    return (sk, ck) if hyp else (sk, - ck)

def t_tan_sec2(t, s2, u, k, hyp=False):
    if k == 0:
        return (tanh(u[0]), 1.0 - tanh(u[0])**2) if hyp else (tan(u[0]), 1.0 + tan(u[0])**2)
    tk = sum((k - j) * s2[j] * u[k - j] for j in range(k)) / k
    sk = 2.0 * (t[0] * tk + sum((k - j) * t[j] * t[k - j] for j in range(1, k)) / k)
    return (tk, - sk) if hyp else (tk, sk)


class Series:

    def __init__(self, jet):
        self.n = len(jet)
        self.jet = [jet[k] for k in range(self.n)]

    @classmethod
    def get(cls, order, value=0.0):
        return cls(t_jet(order, value))

    def __str__(self):
        return ''.join(f"{self.jet[i]:+.9e} " for i in range(self.n))

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
        assert not isinstance(other, Dual)
        if isinstance(other, type(self)):
            return Series([self.jet[k] + other.jet[k] for k in range(self.n)])
        elif isinstance(other, (int, type(self.val))):
            add_jet = [self.jet[k] for k in range(self.n)]
            add_jet[0] += other
            return Series(add_jet)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __radd__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, (int, type(self.val))):
            return self + other
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __sub__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, type(self)):
            return Series([self.jet[k] - other.jet[k] for k in range(self.n)])
        elif isinstance(other, (int, type(self.val))):
            sub_jet = [self.jet[k] for k in range(self.n)]
            sub_jet[0] -= other
            return Series(sub_jet)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rsub__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, (int, type(self.val))):
            return - self + other
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __mul__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, type(self)):
            return Series([t_prod(self.jet, other.jet, k) for k in range(self.n)])
        elif isinstance(other, (int, type(self.val))):
            return Series([self.jet[k] * other for k in range(self.n)])
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rmul__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, (int, type(self.val))):
            return self * other
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __truediv__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, type(self)):
            assert abs(other.val) != 0.0, f"other.val = {other.val}"
            div_jet = t_jet(self.n)
            for k in range(self.n):
                div_jet[k] = t_quot(div_jet, self.jet, other.jet, k)
            return Series(div_jet)
        elif isinstance(other, (int, type(self.val))):
            assert abs(other) != 0.0, f"other = {other}"
            return Series([self.jet[k] / other for k in range(self.n)])
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rtruediv__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, (int, type(self.val))):
            assert abs(self.val) != 0.0, f"self.val = {self.val}"
            return Series(t_jet(self.n, other)) / self
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __pow__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, int):
            i_pow = self
            for _ in range(abs(other) - 1):
                i_pow = i_pow * self
            return i_pow if other > 0 else (Series(t_jet(self.n, 1.0)) / i_pow if other < 0 else Series(t_jet(self.n, 1.0)))
        else:
            assert self.val > 0.0, f"self.val = {self.val}"
            if isinstance(other, type(self)):
                return (self.ln * other).exp
            elif isinstance(other, type(self.val)):
                pow_jet = t_jet(self.n)
                for k in range(self.n):
                    pow_jet[k] = t_pwr(pow_jet, self.jet, other, k)
                return Series(pow_jet)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rpow__(self, other):
        assert not isinstance(other, Dual)
        if isinstance(other, (int, type(self.val))):
            assert other > 0.0, f"other = {other}"
            return (log(other) * self).exp
        raise RuntimeError(f"Incompatible Type: {type(other)}")

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
    def tan(self):
        return self._double(t_tan_sec2)[0]

    @property
    def sinh(self):
        return self._double(t_sin_cos, hyp=True)[0]

    @property
    def cosh(self):
        return self._double(t_sin_cos, hyp=True)[1]

    @property
    def tanh(self):
        return self._double(t_tan_sec2, hyp=True)[0]

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
        return cls(value if isinstance(value, float) else float(value), 0.0)

    def __str__(self):
        return f"{self.val:+.9e} {self.der:+.9e}"

    def __abs__(self):
        return Dual(abs(self.val), self.der if self.val > 0.0 else (- self.der if self.val < 0.0 else 0.0))

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, type(self)):
            return Dual(self.val + other.val, self.der + other.der)
        elif isinstance(other, (int, type(self.val))):
            return Dual(self.val + other, self.der)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __radd__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, (int, type(self.val))):
            return self + other
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __sub__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, type(self)):
            return Dual(self.val - other.val, self.der - other.der)
        elif isinstance(other, (int, type(self.val))):
            return Dual(self.val - other, self.der)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rsub__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, (int, type(self.val))):
            return Dual(- self.val + other, - self.der)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __mul__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, type(self)):
            return Dual(self.val * other.val, self.der * other.val + self.val * other.der)
        elif isinstance(other, (int, type(self.val))):
            return Dual(self.val * other, self.der * other)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rmul__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, (int, type(self.val))):
            return self * other
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __truediv__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, type(self)):
            assert abs(other.val) != 0.0, f"other.val = {other.val}"
            return Dual(self.val / other.val, (self.der * other.val - self.val * other.der) / other.val**2)
        elif isinstance(other, (int, type(self.val))):
            assert abs(other) != 0.0, f"other = {other}"
            return Dual(self.val / other, self.der / other)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rtruediv__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, (int, type(self.val))):
            assert abs(self.val) != 0.0, f"self.val = {self.val}"
            return Dual(other / self.val, - other * self.der / self.val**2)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __pow__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, int):
            i_pow = self
            for _ in range(abs(other) - 1):
                i_pow = i_pow * self
            return i_pow if other > 0 else (Dual(1.0, 0.0) / i_pow if other < 0 else Dual(1.0, 0.0))
        else:
            assert self.val > 0.0, f"self.val = {self.val}"
            if isinstance(other, type(self)):
                return (self.ln * other).exp
            elif isinstance(other, type(self.val)):
                return Dual(self.val**other, other * self.val**(other - 1) * self.der)
        raise RuntimeError(f"Incompatible Type: {type(other)}")

    def __rpow__(self, other):
        assert not isinstance(other, Series)
        if isinstance(other, (int, type(self.val))):
            assert other > 0.0, f"other = {other}"
            return (log(other) * self).exp
        raise RuntimeError(f"Incompatible Type: {type(other)}")

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
    def var(self):
        return Dual(self.val, 1.0)


print(__name__ + " module loaded", file=stderr)
