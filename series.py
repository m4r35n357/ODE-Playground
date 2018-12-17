#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from taylor import t_jet, t_prod, t_quot, t_sqr, t_exp, t_sin_cos, t_tan_sec2, t_pwr, t_ln, t_sqrt, t_asin, t_acos, \
    t_atan


class Series:

    def __init__(self, jet, diff=False):
        self.jet = jet
        self.n = len(self.jet)
        if diff:
            self.jet[1] = 1.0

    def __str__(self):
        string = ""
        for i in range(0, self.n):
            string += "{:.6e} ".format(self.jet[i], end='')
        return string

    def __abs__(self):
        return self.sqr.sqrt

    def __pos__(self):
        return Series([self.jet[k] for k in range(self.n)])

    def __neg__(self):
        return Series([- self.jet[k] for k in range(self.n)])

    def __invert__(self):  # overload - returns a derivative Series
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
            assert not isinstance(other, complex)
            add_jet = [self.jet[k] for k in range(self.n)]
            add_jet[0] += other
            return Series(add_jet)

    def __radd__(self, other):
        assert not isinstance(other, complex)
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            return Series([self.jet[k] - other.jet[k] for k in range(self.n)])
        else:
            assert not isinstance(other, complex)
            sub_jet = [self.jet[k] for k in range(self.n)]
            sub_jet[0] -= other
            return Series(sub_jet)

    def __rsub__(self, other):
        assert not isinstance(other, complex)
        sub_jet = [- self.jet[k] for k in range(self.n)]
        sub_jet[0] += other
        return Series(sub_jet)

    def __mul__(self, other):
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            return Series([t_prod(self.jet, other.jet, k) for k in range(self.n)])
        else:
            assert not isinstance(other, complex)
            return Series([self.jet[k] * other for k in range(self.n)])

    def __rmul__(self, other):
        assert not isinstance(other, complex)
        return self.__mul__(other)

    def __truediv__(self, other):
        div_jet = t_jet(self.n)
        if isinstance(other, Series):
            assert len(other.jet) == self.n
            for k in range(self.n):
                div_jet[k] = t_quot(div_jet, self.jet, other.jet, k)
        else:
            assert not isinstance(other, complex)
            for k in range(self.n):
                div_jet[k] = self.jet[k] / other
        return Series(div_jet)

    def __rtruediv__(self, other):
        assert not isinstance(other, complex)
        other_jet = t_jet(self.n, other)
        rdiv_jet = t_jet(self.n)
        for k in range(self.n):
            rdiv_jet[k] = t_quot(rdiv_jet, other_jet, self.jet, k)
        return Series(rdiv_jet)

    def __pow__(self, a):
        assert not isinstance(a, complex)
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
        return Series([t_sqr(self.jet, k) for k in range(self.n)])

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


print(__name__ + " module loaded", file=stderr)
