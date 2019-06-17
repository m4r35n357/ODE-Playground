#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from sys import stderr
from gmpy2 import mpfr, sqrt, sin_cos, sin, cos, sinh_cosh, sinh, cosh, tan, sec, tanh, atan, asin, acos, exp, log

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

t_ln = lambda l, u, k: log(u[0]) if k == 0 else (u[k] - ddot(u, l, k)) / u[0]

t_pwr = lambda p, u, a, k: u[0]**a if k == 0 else (a * (p[0] * u[k] + ddot(p, u, k)) - ddot(u, p, k)) / u[0]

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
        sec2 = 2 * (t[0] * tn + ddot(t, t, k))
        return tn, - sec2 if hyp else sec2

def t_atan(h, v, u, k):
    if k == 0:
        return atan(u[0]), 1 + u[0]**2
    else:
        return (u[k] - ddot(v, h, k)) / v[0], 2 * (u[0] * u[k] + ddot(u, u, k))

def t_asin(h, v, u, k):
    if k == 0:
        return asin(u[0]), cos(asin(u[0]))
    else:
        hk = (u[k] - ddot(v, h, k)) / v[0]
        return hk, - u[0] * hk - ddot(u, h, k)

def t_acos(h, v, u, k):
    if k == 0:
        return acos(u[0]), - sin(acos(u[0]))
    else:
        hk = (u[k] + ddot(v, h, k)) / v[0]
        return hk, u[0] * hk + ddot(u, h, k)


class Series:

    def __init__(self, jet, variable=False):
        self._jet = jet
        self._n = len(self.jet)
        if variable:
            self._jet[1] = to_mpfr(1)

    @classmethod
    def get(cls, order, value=0, variable=False):
        assert isinstance(order, int)
        return cls(t_jet(order, value), variable)

    def __str__(self):
        string = ""
        for i in range(0, self.n):
            string += f"{self.jet[i]:+.6e} "
        return string

    def __abs__(self):
        return self.sqr.sqrt

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
            assert len(other.jet) == self.n
            return Series([self.jet[k] + other.jet[k] for k in range(self.n)])
        else:
            add_jet = [self.jet[k] for k in range(self.n)]
            add_jet[0] += other
            return Series(add_jet)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            return Series([self.jet[k] - other.jet[k] for k in range(self.n)])
        else:
            sub_jet = [self.jet[k] for k in range(self.n)]
            sub_jet[0] -= other
            return Series(sub_jet)

    def __rsub__(self, other):
        sub_jet = [- self.jet[k] for k in range(self.n)]
        sub_jet[0] += other
        return Series(sub_jet)

    def __mul__(self, other):
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            return Series([t_prod(self.jet, other.jet, k) for k in range(self.n)])
        else:
            return Series([self.jet[k] * other for k in range(self.n)])

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        div_jet = t_jet(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                div_jet[k] = t_quot(div_jet, self.jet, other.jet, k)
        else:
            for k in range(self.n):
                div_jet[k] = self.jet[k] / other
        return Series(div_jet)

    def __rtruediv__(self, other):
        other_jet = t_jet(self.n, other)
        rdiv_jet = t_jet(self.n)
        for k in range(self.n):
            rdiv_jet[k] = t_quot(rdiv_jet, other_jet, self.jet, k)
        return Series(rdiv_jet)

    def __pow__(self, a):
        pow_jet = t_jet(self.n)
        for k in range(self.n):
            pow_jet[k] = t_pwr(pow_jet, self.jet, a, k)
        return Series(pow_jet)

    def _trans(self, fun):
        jet = t_jet(self.n)
        for k in range(self.n):
            jet[k] = fun(jet, self.jet, k)
        return Series(jet)

    def _geom(self, fun, hyp):
        jet_a, jet_b = t_jet(self.n), t_jet(self.n)
        for k in range(self.n):
            jet_a[k], jet_b[k] = fun(jet_a, jet_b, self.jet, k, hyp)
        return Series(jet_a), Series(jet_b)

    def _arc(self, fun):
        h_jet, v_jet = t_jet(self.n), t_jet(self.n)
        for k in range(self.n):
            h_jet[k], v_jet[k] = fun(h_jet, v_jet, self.jet, k)
        return Series(h_jet)

    @property
    def sqr(self):
        return self * self

    @property
    def sqrt(self):
        return self._trans(t_sqrt)

    @property
    def exp(self):
        return self._trans(t_exp)

    @property
    def ln(self):
        return self._trans(t_ln)

    @property
    def sin(self):
        return self._geom(t_sin_cos, False)[0]

    @property
    def cos(self):
        return self._geom(t_sin_cos, False)[1]

    @property
    def sin_cos(self):
        return self._geom(t_sin_cos, False)

    @property
    def sinh(self):
        return self._geom(t_sin_cos, True)[0]

    @property
    def cosh(self):
        return self._geom(t_sin_cos, True)[1]

    @property
    def sinh_cosh(self):
        return self._geom(t_sin_cos, True)

    @property
    def tan(self):
        return self._geom(t_tan_sec2, False)[0]

    @property
    def tanh(self):
        return self._geom(t_tan_sec2, True)[0]

    @property
    def asin(self):
        return self._arc(t_asin)

    @property
    def acos(self):
        return self._arc(t_acos)

    @property
    def atan(self):
        return self._arc(t_atan)

    @property
    def val(self):
        return self.jet[0]

    @val.setter
    def val(self, value):
        self.jet[0] = value

    @property
    def n(self):
        return self._n

    @property
    def jet(self):
        return self._jet


class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value=0, variable=False):
        return cls(to_mpfr(value), to_mpfr(1) if variable else to_mpfr(0))

    def __str__(self):
        return f"{self.val:+.6e} {self.der:+.6e}"

    def __abs__(self):
        return self.sqr.sqrt

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val + other.val, self.der + other.der)
        else:
            return Dual(self.val + other, self.der)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val - other.val, self.der - other.der)
        else:
            return Dual(self.val - other, self.der)

    def __rsub__(self, other):
        return Dual(- self.val + other, - self.der)

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val * other.val, self.der * other.val + self.val * other.der)
        else:
            return Dual(self.val * other, self.der * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val / other.val, (self.der * other.val - self.val * other.der) / other.val**2)
        else:
            return Dual(self.val / other, self.der / other)

    def __rtruediv__(self, other):
        return Dual(other / self.val, - other * self.der / self.val**2)

    def __pow__(self, a):
        return Dual(self.val**a, a * self.val**(a - 1) * self.der)

    @property
    def sqr(self):
        return Dual(self.val**2, 2 * self.der * self.val)

    @property
    def sqrt(self):
        sqrt_val = sqrt(self.val)
        return Dual(sqrt_val, self.der / (2 * sqrt_val))

    @property
    def exp(self):
        exp_val = exp(self.val)
        return Dual(exp_val, self.der * exp_val)

    @property
    def ln(self):
        return Dual(log(self.val), self.der / self.val)

    @property
    def sin(self):
        return Dual(sin(self.val), self.der * cos(self.val))

    @property
    def cos(self):
        return Dual(cos(self.val), - self.der * sin(self.val))

    @property
    def sin_cos(self):
        s, c = sin_cos(self.val)
        return Dual(s, self.der * c), Dual(c, - self.der * s)

    @property
    def sinh(self):
        return Dual(sinh(self.val), self.der * cosh(self.val))

    @property
    def cosh(self):
        return Dual(cosh(self.val), self.der * sinh(self.val))

    @property
    def sinh_cosh(self):
        s, c = sinh_cosh(self.val)
        return Dual(s, self.der * c), Dual(c, self.der * s)

    @property
    def tan(self):
        return Dual(tan(self.val), self.der * sec(self.val)**2)

    @property
    def asin(self):
        return Dual(asin(self.val), self.der / sqrt(1 - self.val**2))

    @property
    def acos(self):
        return Dual(acos(self.val), - self.der / sqrt(1 - self.val**2))

    @property
    def atan(self):
        return Dual(atan(self.val), self.der / (1 + self.val**2))

    @property
    def var(self):
        return Dual(self.val, to_mpfr(1))


print(__name__ + " module loaded", file=stderr)
