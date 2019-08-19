#!/usr/bin/env python3
#  Mutation Testing
#  mut.py --runner pytest --target ad.py --unit-test ad_test -c --disable-operator AOR CRP DDL CDI SDI SDL SVD -e
#  mutmut --paths-to-mutate ad.py run --runner 'pytest ad_test.py'; mutmut results
#  cosmic-ray init config.toml my_session.sqlite
#  cosmic-ray exec my_session.sqlite
#  cr-html my_session.sqlite > my_session.html
from math import pi, e, exp, log, sin, cos, tan, sinh, cosh, tanh
from ad import t_jet, t_horner, Series, Dual

order = 6
ε = 1e-12
zero = 0.0
a = 3.0
b = 4.0
c = 0.5
d = 5

zero_d, zero_s = Dual.get().var, Series.get(order).var
u, v = Dual.get(c).var, Series.get(order, c).var
w, x = Dual.get(b).var, Series.get(order, b).var
y, z = Dual.get(a).var, Series.get(order, a).var

def test_t_jet():
    jet = t_jet(order)
    for k in range(z.n):
        assert abs(jet[k]) < ε
    jet = t_jet(order, c)
    assert abs(jet[0] - c) < ε
    for k in range(1, z.n):
        assert abs(jet[k]) < ε
    jet = t_jet(order, - c)
    assert abs(jet[0] + c) < ε
    for k in range(1, z.n):
        assert abs(jet[k]) < ε

def test_horner():
    assert abs(t_horner([-19, 7, -4, 6], 3) - 128) < ε
    assert abs(t_horner([-19.0, 7.0, -4.0, 6.0], 3.0) - 128.0) < ε

def test_unary_plus():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = + y
    series = ~(+ z)
    assert abs(dual.val - a) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - a) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_unary_minus():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = - y
    series = ~(- z)
    assert abs(dual.val + a) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val + a) < ε
    assert abs(series.jet[1] + 1.0) < ε

def test_abs():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = abs(u)
    series = ~(abs(v))
    assert abs(dual.val - c) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - c) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for k in range(2, z.n):
        assert abs(series.jet[k]) < ε
    dual = abs(- u)
    series = ~(abs(- v))
    assert abs(dual.val - c) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - c) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for k in range(2, z.n):
        assert abs(series.jet[k]) < ε
    dual = abs(zero_d)
    series = ~(abs(zero_s))
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_add_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y + w
    series = ~(z + x)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 2.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 2.0) < ε

def test_add_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y + d
    series = ~(z + d)
    assert abs(dual.val - (a + d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + d)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_add_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d + y
    series = ~(d + z)
    assert abs(dual.val - (a + d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + d)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_add_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y + b
    series = ~(z + b)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_add_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b + y
    series = ~(b + z)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_subtract_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y - w
    series = ~(z - x)
    assert abs(dual.val - (a - b)) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - (a - b)) < ε
    assert abs(series.jet[1]) < ε

def test_subtract_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y - d
    series = ~(z - d)
    assert abs(dual.val - (a - d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a - d)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_subtract_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d - y
    series = ~(d - z)
    assert abs(dual.val - (d - a)) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val - (d - a)) < ε
    assert abs(series.jet[1] + 1.0) < ε

def test_subtract_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y - b
    series = ~(z - b)
    assert abs(dual.val - (a - b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a - b)) < ε
    assert abs(series.jet[1] - 1.0) < ε

def test_subtract_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b - y
    series = ~(b - z)
    assert abs(dual.val - (b - a)) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val - (b - a)) < ε
    assert abs(series.jet[1] + 1.0) < ε

def test_multiply_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y * w
    series = ~(z * x)
    assert abs(dual.val - a * b) < ε
    assert abs(dual.der - (a + b)) < ε
    assert abs(series.val - a * b) < ε
    assert abs(series.jet[1] - (a + b)) < ε

def test_multiply_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y * d
    series = ~(z * d)
    assert abs(dual.val - a * d) < ε
    assert abs(dual.der - d) < ε
    assert abs(series.val - a * d) < ε
    assert abs(series.jet[1] - d) < ε

def test_multiply_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d * y
    series = ~(d * z)
    assert abs(dual.val - d * a) < ε
    assert abs(dual.der - d) < ε
    assert abs(series.val - d * a) < ε
    assert abs(series.jet[1] - d) < ε

def test_multiply_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y * b
    series = ~(z * b)
    assert abs(dual.val - a * b) < ε
    assert abs(dual.der - b) < ε
    assert abs(series.val - a * b) < ε
    assert abs(series.jet[1] - b) < ε

def test_multiply_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b * y
    series = ~(b * z)
    assert abs(dual.val - b * a) < ε
    assert abs(dual.der - b) < ε
    assert abs(series.val - b * a) < ε
    assert abs(series.jet[1] - b) < ε

def test_divide_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y / w
    series = ~(z / x)
    assert abs(dual.val - a / b) < ε
    assert abs(dual.der - (b - a) / b**2) < ε
    assert abs(series.val - a / b) < ε
    assert abs(series.jet[1] - (b - a) / b**2) < ε

def test_divide_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y / d
    series = ~(z / d)
    assert abs(dual.val - a / d) < ε
    assert abs(dual.der - 1.0 / d) < ε
    assert abs(series.val - a / d) < ε
    assert abs(series.jet[1] - 1.0 / d) < ε

def test_divide_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d / y
    series = ~(d / z)
    assert abs(dual.val - d / a) < ε
    assert abs(dual.der + d / a**2) < ε
    assert abs(series.val - d / a) < ε
    assert abs(series.jet[1] + d / a**2) < ε

def test_divide_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y / b
    series = ~(z / b)
    assert abs(dual.val - a / b) < ε
    assert abs(dual.der - 1.0 / b) < ε
    assert abs(series.val - a / b) < ε
    assert abs(series.jet[1] - 1.0 / b) < ε

def test_divide_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b / y
    series = ~(b / z)
    assert abs(dual.val - b / a) < ε
    assert abs(dual.der + b / a**2) < ε
    assert abs(series.val - b / a) < ε
    assert abs(series.jet[1] + b / a**2) < ε

def test_reciprocal():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 1.0 / y
    series = ~(1.0 / z)
    assert abs(dual.val - 1.0 / a) < ε
    assert abs(dual.der + 1.0 / a**2) < ε
    assert abs(series.val - 1.0 / a) < ε
    assert abs(series.jet[1] + 1.0 / a**2) < ε

def test_pow_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y**w
    series = ~(z**x)
    assert abs(dual.val - a**b) < ε
    assert abs(series.val - a**b) < ε

def test_pow_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d**y
    series = ~d**z
    assert abs(dual.val - d**a) < ε
    assert abs(series.val - d**a) < ε

def test_pow_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b**y
    series = ~b**z
    assert abs(dual.val - b**a) < ε
    assert abs(series.val - b**a) < ε

def test_pow_int():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(d, int)
    dual = w**d
    series = ~x**d
    assert abs(dual.val - b**d) < ε
    assert abs(dual.der - d * b**(d - 1)) < ε
    assert abs(series.val - b**d) < ε
    assert abs(series.jet[1] - d * b**(d - 1)) < ε

def test_pow_int_neg():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(d, int)
    dual = w**-d
    series = ~x**-d
    assert abs(dual.val - 1.0 / b**d) < ε
    assert abs(series.val - 1.0 / b**d) < ε

def test_pow_float():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(a, float)
    dual = w**a
    series = ~x**a
    assert abs(dual.val - b**a) < ε
    assert abs(dual.der - a * b**(a - 1.0)) < ε
    assert abs(series.val - b**a) < ε
    assert abs(series.jet[1] - a * b**(a - 1.0)) < ε

def test_pow_float_neg():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(a, float)
    dual = w**-a
    series = ~x**-a
    assert abs(dual.val - 1.0 / b**a) < ε
    assert abs(series.val - 1.0 / b**a) < ε

def test_pow_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w**0
    series = ~x**0
    assert abs(dual.val - 1.0) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - 1.0) < ε
    assert abs(series.jet[1]) < ε

def test_pow_zero_float():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w**0.0
    series = ~x**0.0
    assert abs(dual.val - 1.0) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - 1.0) < ε
    assert abs(series.jet[1]) < ε

def test_exp():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y.exp
    series = ~z.exp
    assert abs(dual.val - exp(a)) < ε
    assert abs(dual.der - exp(a)) < ε
    for k in range(z.n):
        assert abs(series.jet[k] - exp(a)) < ε

def test_ln():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y.ln
    series = z.ln
    assert abs(dual.val - log(a)) < ε
    assert abs(dual.der - 1.0 / a) < ε
    assert abs(series.val - log(a)) < ε
    assert abs(series.jet[1] - 1.0 / a) < ε

def test_sin():
    dual = Dual.get(pi / a).var.sin
    series = ~Series.get(order, pi / a).var.sin
    assert abs(dual.val - sin(pi / a)) < ε
    assert abs(dual.der - cos(pi / a)) < ε
    for k in range(z.n):
        if k % 4 == 0:
            assert abs(series.jet[k] - sin(pi / a)) < ε
        elif k % 4 == 1:
            assert abs(series.jet[k] - cos(pi / a)) < ε
        elif k % 4 == 2:
            assert abs(series.jet[k] + sin(pi / a)) < ε
        elif k % 4 == 3:
            assert abs(series.jet[k] + cos(pi / a)) < ε

def test_cos():
    dual = Dual.get(pi / a).var.cos
    series = ~Series.get(order, pi / a).var.cos
    assert abs(dual.val - cos(pi / a)) < ε
    assert abs(dual.der + sin(pi / a)) < ε
    for k in range(z.n):
        if k % 4 == 0:
            assert abs(series.jet[k] - cos(pi / a)) < ε
        elif k % 4 == 1:
            assert abs(series.jet[k] + sin(pi / a)) < ε
        elif k % 4 == 2:
            assert abs(series.jet[k] + cos(pi / a)) < ε
        elif k % 4 == 3:
            assert abs(series.jet[k] - sin(pi / a)) < ε

def test_tan():
    dual = Dual.get(pi / b).var.tan
    series = ~Series.get(order, pi / b).var.tan
    assert abs(dual.val - tan(pi / b)) < ε
    assert abs(dual.der - (1.0 + tan(pi / b)**2)) < ε
    assert abs(series.val - tan(pi / b)) < ε
    assert abs(series.jet[1] - (1.0 + tan(pi / b)**2)) < ε

def test_sinh():
    dual = Dual.get(pi / a).var.sinh
    series = ~Series.get(order, pi / a).var.sinh
    assert abs(dual.val - sinh(pi / a)) < ε
    assert abs(dual.der - cosh(pi / a)) < ε
    assert abs(series.val - sinh(pi / a)) < ε
    assert abs(series.jet[1] - cosh(pi / a)) < ε

def test_cosh():
    dual = Dual.get(pi / a).var.cosh
    series = ~Series.get(order, pi / a).var.cosh
    assert abs(dual.val - cosh(pi / a)) < ε
    assert abs(dual.der - sinh(pi / a)) < ε
    assert abs(series.val - cosh(pi / a)) < ε
    assert abs(series.jet[1] - sinh(pi / a)) < ε

def test_tanh():
    dual = Dual.get(pi / b).var.tanh
    series = ~Series.get(order, pi / b).var.tanh
    assert abs(dual.val - tanh(pi / b)) < ε
    assert abs(dual.der - (1.0 - tanh(pi / b)**2)) < ε
    assert abs(series.val - tanh(pi / b)) < ε
    assert abs(series.jet[1] - (1.0 - tanh(pi / b)**2)) < ε

#  Zero identities

def test_pow_neg_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = 1.0 / w**2 - w**-2.0
    series = ~(1.0 / x**2 - x**-2.0)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_pow_pos_frac_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = (w**2)**0.5 - abs(w)
    series = ~((x**2)**0.5) - abs(x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_exp_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.exp - e**w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_ln_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.exp.ln - w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_sinh_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 0.5 * (y.exp - (-y).exp) - y.sinh
    series = ~(0.5 * (z.exp - (-z).exp) - z.sinh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_cosh_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 0.5 * (y.exp + (-y).exp) - y.cosh
    series = ~(0.5 * (z.exp + (-z).exp) - z.cosh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_tan_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.tan - w.sin / w.cos
    series = ~(x.tan - x.sin / x.cos)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_tanh_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = u.tanh - u.sinh / u.cosh
    series = ~(v.tanh - v.sinh / v.cosh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_sin_3x_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = (3 * y).sin - 3.0 * y.sin + 4.0 * y.sin**3
    series = ~((3 * z).sin - 3.0 * z.sin + 4.0 * z.sin**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_cos_3x_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = (3 * y).cos + 3.0 * y.cos - 4.0 * y.cos**3
    series = ~((3 * z).cos + 3.0 * z.cos - 4.0 * z.cos**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_sinh_3x_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = (3 * u).sinh - 3.0 * u.sinh - 4.0 * u.sinh**3
    series = ~((3 * v).sinh - 3.0 * v.sinh - 4.0 * v.sinh**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε

def test_cosh_3x_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = (3 * u).cosh + 3.0 * u.cosh - 4.0 * u.cosh**3
    series = ~((3 * v).cosh + 3.0 * v.cosh - 4.0 * v.cosh**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for k in range(z.n):
        assert abs(series.jet[k]) < ε
