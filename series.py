
from taylor import jet_0, jet_c, t_prod, t_quot, t_sqr, t_exp, t_sin_cos, t_tan_sec2, t_pwr, t_ln


class Series:

    def __init__(self, jet, diff=False):
        self.jet = jet
        self.n = len(self.jet)
        if diff:
            self.jet[1] = 1.0
        self.diff_status = diff

    def __str__(self):
        string = "{:14.6e} ".format(self.jet[0], end='')
        for i in range(1, self.n):
            string += "{:14.6e}".format(self.jet[i], end='')
        return string

    def __neg__(self):
        jet = jet_0(self.n)
        for k in range(self.n):
            jet[k] = - self.jet[k]
        return Series(jet)

    def __add__(self, other):
        jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                jet[k] = self.jet[k] + other.jet[k]
        else:
            jet[0] = self.jet[0] + other
            for k in range(1, self.n):
                jet[k] = self.jet[k]
        return Series(jet)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                jet[k] = self.jet[k] - other.jet[k]
        else:
            jet[0] = self.jet[0] - other
            for k in range(1, self.n):
                jet[k] = self.jet[k]
        return Series(jet)

    def __rsub__(self, other):
        jet = jet_0(self.n)
        jet[0] = other - self.jet[0]
        for k in range(1, self.n):
            jet[k] = - self.jet[k]
        return Series(jet)

    def __mul__(self, other):
        jet = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                jet[k] = t_prod(self.jet, other.jet, k)
        else:
            for k in range(self.n):
                jet[k] = self.jet[k] * other
        return Series(jet)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        q = jet_0(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                q[k] = t_quot(q, self.jet, other.jet, k)
        else:
            for k in range(self.n):
                q[k] = self.jet[k] / other
        return Series(q)

    def __rtruediv__(self, other):
        q = (self ** -1).jet
        for k in range(self.n):
            q[k] = other * q[k]
        return Series(q)

    def __pow__(self, a):
        pwr_jet = jet_0(self.n)
        for k in range(self.n):
            pwr_jet[k] = t_pwr(pwr_jet, self.jet, a, k)
        return Series(pwr_jet)

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
        return self._t_s2(False)[0]

    @property
    def tan_sec2(self):
        return self._t_s2(False)

    @property
    def tanh(self):
        return self._t_s2(True)[0]

    @property
    def tanh_sech2(self):
        return self._t_s2(True)

    def _t_s2(self, hyp):
        tan_jet = jet_0(self.n)
        sec2_jet = jet_0(self.n)
        for k in range(self.n):
            tan_jet[k], sec2_jet[k] = t_tan_sec2(tan_jet, sec2_jet, self.jet, k, hyp)
        return Series(tan_jet), Series(sec2_jet)

    @property
    def ln(self):
        ln_jet = jet_0(self.n)
        for k in range(self.n):
            ln_jet[k] = t_ln(ln_jet, self.jet, k)
        return Series(ln_jet)


def newton(model, initial, target=0.0, tol=1.0e-12, max_it=100):
    x = Series(jet_c(initial, 2), diff=True)
    t = Series(jet_c(target, 2))
    f = Series(jet_c(1.0, 2))
    delta = 1.0
    counter = 1
    while abs(f.jet[0]) > tol or abs(delta) > tol:
        f = model(x, t)
        delta = - f.jet[0] / f.jet[1]
        x.jet[0] += delta
        print("{:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, x.jet[0], f.jet[0] + t.jet[0], delta))
        counter += 1
        if counter > max_it:
            break
    return counter, x.jet[0], f.jet[0] + t.jet[0], delta


def householder(model, initial, n, target=0.0, tol=1.0e-12, max_it=100):
    x = Series(jet_c(initial, n), diff=True)
    t = Series(jet_c(target, n))
    f = Series(jet_c(1.0, n))
    delta = 1.0
    counter = 1
    while abs(f.jet[0]) > tol or abs(delta) > tol:
        f = model(x, t)
        r = (f ** -1).derivatives
        delta = r.jet[n - 2] / r.jet[n - 1]
        x.jet[0] += delta * (n - 1)
        print("{:3d} {:22.15e} {:22.15e} {:10.3e}".format(counter, x.jet[0], f.jet[0] + t.jet[0], delta))
        counter += 1
        if counter > max_it:
            break
    return counter, x.jet[0], f.jet[0] + t.jet[0], delta
