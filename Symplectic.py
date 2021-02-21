"""
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from sys import stderr
from math import log, exp, sin, cos, tan, cosh, sinh, tanh


class Context:
    places = 3


class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value=0.0):
        # noinspection PyArgumentList
        return cls(value, 0.0)

    def __str__(self):
        return f'{self.val:+.{Context.places}e} {self.der:+.{Context.places}e}'

    def __abs__(self):
        return Dual(abs(self.val), self.der if self.val > 0.0 else (- self.der if self.val < 0.0 else 0.0))

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val + o.val, self.der + o.der)
        elif isinstance(o, (type(self.val), int)):
            return Dual(self.val + o, self.der)
        raise RuntimeError(f"Incompatible Type: {type(o)}, expecting Dual or {type(self.val)}")

    def __radd__(self, o):
        return self.__add__(o)

    def __sub__(self, o):
        return self.__add__(- o)

    def __rsub__(self, o):
        return (- self).__add__(o)

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.val * o.val, self.der * o.val + self.val * o.der)
        elif isinstance(o, (type(self.val), int)):
            return Dual(self.val * o, self.der * o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rmul__(self, o):
        return self.__mul__(o)

    def __truediv__(self, o):
        if isinstance(o, Dual):
            assert abs(o.val) != 0.0, f"other.val = {o.val}"
            return Dual(self.val / o.val, (self.der * o.val - self.val * o.der) / o.val**2)
        elif isinstance(o, (type(self.val), int)):
            assert abs(o) != 0.0, f"other = {o}"
            return self.__mul__(1.0 / o)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rtruediv__(self, o):
        return Dual(o, 0.0).__truediv__(self)

    def __pow__(self, o):
        if isinstance(o, int):
            i_pow = Dual(self.val, self.der)
            for _ in range(abs(o) - 1):
                i_pow = i_pow.__mul__(self)
            return i_pow if o > 0 else (Dual(1.0, 0.0).__truediv__(i_pow) if o < 0 else Dual(1.0, 0.0))
        else:
            if isinstance(o, Dual):
                return (self.ln.__mul__(o)).exp
            elif isinstance(o, type(self.val)):
                assert self.val > 0.0, f"self.val = {self.val}"  # pragma: no mutate
                return Dual(self.val**o, o * self.val**(o - 1) * self.der)
        raise RuntimeError(f"Incompatible Type: {type(o)}")

    def __rpow__(self, o):
        assert o > 0.0, f"other = {o}"
        return (self.__mul__(log(o))).exp

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


class Symplectic(object):

    def __init__(self, model, h, order, scheme):
        self.model = model
        self.h = h
        if scheme == 'yoshida':
            print("Yoshida composition", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = 2.0
        elif scheme == 'suzuki':
            print("Suzuki composition", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = 4.0
        else:
            raise Exception('>>> Composition scheme must be yoshida or suzuki, was "{found}" <<<'.format(found=scheme))
        if order == 'b1':
            print("1st order (Euler-Cromer)", file=stderr)
            self.method = self.euler_cromer
        elif order == 'b2':
            print("2nd order (Stormer-Verlet))", file=stderr)
            self.method = self.second_order
        elif order == 'b4':
            print("4th order (Composed)", file=stderr)
            self.method = self.fourth_order
        elif order == 'f4':
            print("4th order (Forest-Ruth)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = 2.0
            self.method = self.fourth_order_forest_ruth
        elif order == 'b6':
            print("6th order (Composed)", file=stderr)
            self.method = self.sixth_order
        elif order == 'f6':
            print("6th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = 2.0
            self.method = self.sixth_order_forest_ruth
        elif order == 'b8':
            print("8th order (Composed)", file=stderr)
            self.method = self.eightth_order
        elif order == 'f8':
            print("8th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = 2.0
            self.method = self.eightth_order_forest_ruth
        elif order == 'b10':
            print("10th order (Composed)", file=stderr)
            self.method = self.tenth_order
        elif order == 'f10':
            print("10th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = 2.0
            self.method = self.tenth_order_forest_ruth
        elif order == 's4':
            print("4th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = 4.0
            self.method = self.fourth_order_smith
        elif order == 's6':
            print("6th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = 4.0
            self.method = self.sixth_order_smith
        elif order == 's8':
            print("8th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = 4.0
            self.method = self.eightth_order_smith
        elif order == 's10':
            print("10th order (Smith) (Composed)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = 4.0
            self.method = self.tenth_order_smith
        else:
            raise Exception(
                '>>> Integrator must be b1, b2, [bfs]4, [bfs]6, [bfs]8, or [bfs]10, was "{found}" <<<'.format(
                    found=order))
        self.w1 = 1.0 / (self.scheme_root - self.scheme_root ** (1.0 / 9.0))
        self.x1 = 1.0 / (self.scheme_root - self.scheme_root ** (1.0 / 7.0))
        self.y1 = 1.0 / (self.scheme_root - self.scheme_root ** (1.0 / 5.0))
        self.z1 = 1.0 / (self.scheme_root - self.scheme_root ** (1.0 / 3.0))
        self.w0 = 1.0 - self.scheme_root * self.w1
        self.x0 = 1.0 - self.scheme_root * self.x1
        self.y0 = 1.0 - self.scheme_root * self.y1
        self.z0 = 1.0 - self.scheme_root * self.z1
        self.cd_sv = [
            0.5 * h, h
        ]
        self.cd_f4 = [
            0.5 * h * self.z1, h * self.z1, 0.5 * h * (self.z0 + self.z1), h * self.z0
        ]
        self.cd_s4 = [
            0.5 * h * self.z1, h * self.z1, h * self.z1, h * self.z1, 0.5 * h * (self.z1 + self.z0), h * self.z0
        ]
        self.cd_s6 = [
            0.5 * h * self.z1 * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            0.5 * h * (self.z1 + self.z0) * self.y1, h * self.z0 * self.y1, 0.5 * h * (self.z0 + self.z1) * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            h * self.z1 * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            0.5 * h * (self.z1 + self.z0) * self.y1, h * self.z0 * self.y1, 0.5 * h * (self.z0 + self.z1) * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            0.5 * h * self.z1 * (self.y1 + self.y0),
            h * self.z1 * self.y0, h * self.z1 * self.y0, h * self.z1 * self.y0,
            0.5 * h * (self.z1 + self.z0) * self.y0, h * self.z0 * self.y0
        ]
        self.cd_s8 = [
            0.5 * h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y0 * self.x1, h * self.z0 * self.y0 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y0 * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            0.5 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y0 * self.x1, h * self.z0 * self.y0 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y0 * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            0.5 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            0.5 * h * self.z1 * self.y1 * (self.x1 + self.x0),
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x0, h * self.z0 * self.y1 * self.x0, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            0.5 * h * (self.z1 + self.z0) * self.y1 * self.x0, h * self.z0 * self.y1 * self.x0, 0.5 * h * (self.z0 + self.z1) * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            0.5 * h * self.z1 * (self.y1 + self.y0) * self.x0,
            h * self.z1 * self.y0 * self.x0, h * self.z1 * self.y0 * self.x0, h * self.z1 * self.y0 * self.x0,
            0.5 * h * (self.z1 + self.z0) * self.y0 * self.x0, h * self.z0 * self.y0 * self.x0
        ]

    def euler_cromer(self):
        self.model.q_update(self.h)
        self.model.p_update(self.h)

    def stormer_verlet(self, s):
        self.model.q_update(s * self.cd_sv[0])
        self.model.p_update(s * self.cd_sv[1])
        self.model.q_update(s * self.cd_sv[0])

    def second_order(self):
        self.stormer_verlet(1.0)

    @staticmethod
    def yoshida(base_method, s, plus, minus):
        base_method(s * plus)
        base_method(s * minus)
        base_method(s * plus)

    @staticmethod
    def suzuki(base_method, s, plus, minus):
        base_method(s * plus)
        base_method(s * plus)
        base_method(s * minus)
        base_method(s * plus)
        base_method(s * plus)

    def base4(self, s):
        self.scheme(self.stormer_verlet, s, self.z1, self.z0)

    def base6(self, s):
        self.scheme(self.base4, s, self.y1, self.y0)

    def base8(self, s):
        self.scheme(self.base6, s, self.x1, self.x0)

    def fourth_order(self):
        self.base4(1.0)

    def sixth_order(self):
        self.base6(1.0)

    def eightth_order(self):
        self.base8(1.0)

    def tenth_order(self):
        self.scheme(self.base8, 1.0, self.w1, self.w0)

    def forest_ruth_4(self, s):
        self.model.q_update(s * self.cd_f4[0])
        self.model.p_update(s * self.cd_f4[1])
        self.model.q_update(s * self.cd_f4[2])
        self.model.p_update(s * self.cd_f4[3])
        self.model.q_update(s * self.cd_f4[2])
        self.model.p_update(s * self.cd_f4[1])
        self.model.q_update(s * self.cd_f4[0])

    def fourth_order_forest_ruth(self):
        self.forest_ruth_4(1.0)

    def base6_forest_ruth(self, s):
        self.scheme(self.forest_ruth_4, s, self.y1, self.y0)

    def sixth_order_forest_ruth(self):
        self.base6_forest_ruth(1.0)

    def base8_forest_ruth(self, s):
        self.scheme(self.base6_forest_ruth, s, self.x1, self.x0)

    def eightth_order_forest_ruth(self):
        self.base8_forest_ruth(1.0)

    def tenth_order_forest_ruth(self):
        self.scheme(self.base8_forest_ruth, 1.0, self.w1, self.w0)

    def smith(self, s):
        size = len(self.coefficients)
        for i in range(size):
            (self.model.q_update if i % 2 == 0 else self.model.p_update)(s * self.coefficients[i])
        for i in range(size - 2, -1, -1):
            (self.model.q_update if i % 2 == 0 else self.model.p_update)(s * self.coefficients[i])

    def fourth_order_smith(self):
        self.coefficients = self.cd_s4
        self.smith(1.0)

    def sixth_order_smith(self):
        self.coefficients = self.cd_s6
        self.smith(1.0)

    def eightth_order_smith(self):
        self.coefficients = self.cd_s8
        self.smith(1.0)

    def tenth_order_smith(self):
        self.coefficients = self.cd_s8
        self.scheme(self.smith, 1.0, self.w1, self.w0)


print(__name__ + " module loaded", file=stderr)
