#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from enum import Enum
from taylor import jet_0, jet_c, t_prod, t_quot, t_sqr, t_exp, t_sin_cos, t_tan_sec2, t_pwr, t_ln, t_sqrt, t_asin, \
    t_acos, t_atan


class Solver(Enum):
    ROOT = 0
    EXTREMUM = 1
    INFLECTION = 2


class Series:

    def __init__(self, jet, diff=False):
        self.jet = jet
        self.n = len(self.jet)
        if diff:
            self.jet[1] = 1.0
        self.diff_status = diff

    def __str__(self):
        string = ""
        for i in range(0, self.n):
            string += "{:.6e} ".format(self.jet[i], end='')
        return string

    def __abs__(self):
        return self.sqr.sqrt

    def __neg__(self):
        return Series([- self.jet[k] for k in range(self.n)])

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
        div_jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                div_jet[k] = t_quot(div_jet, self.jet, other.jet, k)
        else:
            for k in range(self.n):
                div_jet[k] = self.jet[k] / other
        return Series(div_jet)

    def __rtruediv__(self, other):
        other_jet = jet_c(other, self.n)
        rdiv_jet = jet_0(self.n)
        for k in range(self.n):
            rdiv_jet[k] = t_quot(rdiv_jet, other_jet, self.jet, k)
        return Series(rdiv_jet)

    def __pow__(self, a):
        pow_jet = jet_0(self.n)
        for k in range(self.n):
            pow_jet[k] = t_pwr(pow_jet, self.jet, a, k)
        return Series(pow_jet)

    def _s_c(self, hyp):
        sin_jet, cos_jet = jet_0(self.n), jet_0(self.n)
        for k in range(self.n):
            sin_jet[k], cos_jet[k] = t_sin_cos(sin_jet, cos_jet, self.jet, k, hyp)
        return Series(sin_jet), Series(cos_jet)

    def _t_s2(self, hyp):
        tan_jet, sec2_jet = jet_0(self.n), jet_0(self.n)
        for k in range(self.n):
            tan_jet[k], sec2_jet[k] = t_tan_sec2(tan_jet, sec2_jet, self.jet, k, hyp)
        return Series(tan_jet)

    def _arc(self, fun):
        h_jet, v_jet = jet_0(self.n), jet_0(self.n)
        for k in range(self.n):
            h_jet[k], v_jet[k] = fun(h_jet, v_jet, self.jet, k)
        return Series(h_jet)

    @property
    def derivatives(self):
        d = jet_c(self.jet[0], self.n)
        fac = 1
        for i in range(1, self.n):
            fac *= i
            d[i] = fac * self.jet[i]
        return Series(d)

    @property
    def sqr(self):
        return Series([t_sqr(self.jet, k) for k in range(self.n)])

    @property
    def sqrt(self):
        sqrt_jet = jet_0(self.n)
        for k in range(self.n):
            sqrt_jet[k] = t_sqrt(sqrt_jet, self.jet, k)
        return Series(sqrt_jet)

    @property
    def exp(self):
        exp_jet = jet_0(self.n)
        for k in range(self.n):
            exp_jet[k] = t_exp(exp_jet, self.jet, k)
        return Series(exp_jet)

    @property
    def sin(self):
        return self._s_c(False)[0]

    @property
    def cos(self):
        return self._s_c(False)[1]

    @property
    def sin_cos(self):
        return self._s_c(False)

    @property
    def sinh(self):
        return self._s_c(True)[0]

    @property
    def cosh(self):
        return self._s_c(True)[1]

    @property
    def sinh_cosh(self):
        return self._s_c(True)

    @property
    def tan(self):
        return self._t_s2(False)

    @property
    def tanh(self):
        return self._t_s2(True)

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
    def ln(self):
        ln_jet = jet_0(self.n)
        for k in range(self.n):
            ln_jet[k] = t_ln(ln_jet, self.jet, k)
        return Series(ln_jet)


def bisect(model, ax, bx, f_tol, x_tol, target=0.0, max_it=100, mode=Solver.ROOT):
    a = Series(jet_c(ax, 3), diff=True)
    b = Series(jet_c(bx, 3), diff=True)
    c = Series(jet_0(3))
    fc = Series(jet_c(1.0, 3))
    f_sign = model(a, target).derivatives
    delta = 1.0
    counter = 1
    while abs(fc.jet[0 + mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / 2.0
        fc = model(c, target).derivatives
        if f_sign.jet[0 + mode.value] * fc.jet[0 + mode.value] < 0.0:
            b = c
        else:
            a = c
        delta = b.jet[0 + mode.value] - a.jet[0 + mode.value]
        counter += 1
        if counter > max_it:
            break
    return counter, c.jet[0], fc.jet[0 + mode.value] + target, delta


def newton(model, initial, f_tol, x_tol, target=0.0, max_it=100, mode=Solver.ROOT):
    x = Series(jet_c(initial, 2 + mode.value), diff=True)
    f = Series(jet_c(1.0, 2 + mode.value))
    delta = 1.0
    counter = 1
    while abs(f.jet[0 + mode.value]) > f_tol or abs(delta) > x_tol:
        f = model(x, target).derivatives
        delta = - f.jet[0 + mode.value] / f.jet[1 + mode.value]
        x.jet[0] += delta
        counter += 1
        if counter > max_it:
            break
    return counter, x.jet[0], f.jet[0 + mode.value] + target, delta


def householder(model, initial, n, f_tol, x_tol, target=0.0, max_it=100, mode=Solver.ROOT):
    x = Series(jet_c(initial, n + mode.value), diff=True)
    f = Series(jet_c(1.0, n + mode.value))
    delta = 1.0
    counter = 1
    while abs(f.jet[0] + mode.value) > f_tol or abs(delta) > x_tol:
        f = model(x, target)
        r = (1 / f).derivatives
        delta = r.jet[n - 2 + mode.value] / r.jet[n - 1 + mode.value]
        x.jet[0] += delta * (n - 1 + mode.value)
        counter += 1
        if counter > max_it:
            break
    return counter, x.jet[0], f.jet[0] + mode.value + target, delta
