#
#  (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr
from gmpy2 import sqrt, sin, cos, mpfr, exp, log, tan, sec, atan, asin, acos


# noinspection PyArgumentList
def make_mpfr(x):
    return mpfr(str(x)) if isinstance(x, (float, int)) else mpfr(x)


class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value, variable=False):
        return cls(make_mpfr(value), make_mpfr(1) if variable else make_mpfr(0))

    def __str__(self):
        return "{:+.6e} {:+.6e}".format(self.val, self.der)

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
    def var(self):
        return Dual(self.val, make_mpfr(1))

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


print(__name__ + " module loaded", file=stderr)
