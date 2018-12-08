from sys import stderr
from enum import Enum
from taylor import jet_0, jet_c, t_prod, t_quot, t_sqr, t_exp, t_sin_cos, t_tan_sec2, t_pwr, t_ln, t_sqrt


class SolveMode(Enum):
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

    def __neg__(self):
        neg_jet = jet_0(self.n)
        for k in range(self.n):
            neg_jet[k] = - self.jet[k]
        return Series(neg_jet)

    def __add__(self, other):
        add_jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                add_jet[k] = self.jet[k] + other.jet[k]
        else:
            add_jet[0] = self.jet[0] + other
            for k in range(1, self.n):
                add_jet[k] = self.jet[k]
        return Series(add_jet)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        sub_jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                sub_jet[k] = self.jet[k] - other.jet[k]
        else:
            sub_jet[0] = self.jet[0] - other
            for k in range(1, self.n):
                sub_jet[k] = self.jet[k]
        return Series(sub_jet)

    def __rsub__(self, other):
        rsub_jet = jet_0(self.n)
        rsub_jet[0] = other - self.jet[0]
        for k in range(1, self.n):
            rsub_jet[k] = - self.jet[k]
        return Series(rsub_jet)

    def __mul__(self, other):
        mul_jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                mul_jet[k] = t_prod(self.jet, other.jet, k)
        else:
            for k in range(self.n):
                mul_jet[k] = self.jet[k] * other
        return Series(mul_jet)

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

    @property
    def diff(self):
        return self.diff_status

    @diff.setter
    def diff(self, value):
        self.diff_status = value

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
        sqr_jet = jet_0(self.n)
        for k in range(self.n):
            sqr_jet[k] = t_sqr(self.jet, k)
        return Series(sqr_jet)

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

    def _s_c(self, hyp):
        sin_jet = jet_0(self.n)
        cos_jet = jet_0(self.n)
        for k in range(self.n):
            sin_jet[k], cos_jet[k] = t_sin_cos(sin_jet, cos_jet, self.jet, k, hyp)
        return Series(sin_jet), Series(cos_jet)

    @property
    def tan(self):
        return self._t_s2(False)

    @property
    def tanh(self):
        return self._t_s2(True)

    def _t_s2(self, hyp):
        tan_jet = jet_0(self.n)
        sec2_jet = jet_0(self.n)
        for k in range(self.n):
            tan_jet[k], sec2_jet[k] = t_tan_sec2(tan_jet, sec2_jet, self.jet, k, hyp)
        return Series(tan_jet)

    @property
    def ln(self):
        ln_jet = jet_0(self.n)
        for k in range(self.n):
            ln_jet[k] = t_ln(ln_jet, self.jet, k)
        return Series(ln_jet)

    @property
    def abs(self):
        return self.sqr.sqrt


def bisect(model, ax, bx, f_tol, x_tol, target=0.0, max_it=100, mode=SolveMode.ROOT):
    a = Series(jet_c(ax, 1 + mode.value))
    b = Series(jet_c(bx, 1 + mode.value))
    c = Series(jet_0(1 + mode.value))
    fa = model(a, target)
    fc = Series(jet_c(1.0, 1 + mode.value))
    delta = 1.0
    counter = 1
    while abs(fc.jet[0 + mode.value]) > f_tol or abs(delta) > x_tol:
        c = (a + b) / 2.0
        fc = model(c, target).derivatives
        if fa.jet[0 + mode.value] * fc.jet[0 + mode.value] < 0.0:
            b = c
        else:
            a = c
        delta = b.jet[0 + mode.value] - a.jet[0 + mode.value]
        counter += 1
        if counter > max_it:
            break
    print("    {:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, c.jet[0], fc.jet[0 + mode.value] + target, delta),
          file=stderr)
    return counter, c.jet[0], fc.jet[0 + mode.value] + target, delta


def newton(model, initial, f_tol, x_tol, target=0.0, max_it=100, mode=SolveMode.ROOT):
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
    print("    {:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, x.jet[0], f.jet[0 + mode.value] + target, delta),
          file=stderr)
    return counter, x.jet[0], f.jet[0 + mode.value] + target, delta


def householder(model, initial, n, f_tol, x_tol, target=0.0, max_it=100, mode=SolveMode.ROOT):
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
    print("    {:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, x.jet[0], f.jet[0] + mode.value + target, delta),
          file=stderr)
    return counter, x.jet[0], f.jet[0] + mode.value + target, delta
