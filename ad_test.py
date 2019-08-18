#!/usr/bin/env python3

from math import pi, e
from ad import Series, Dual

order = 6
ε=1e-15
a = 3.0
b = 4.0
c = 0.5
d = 3

u, v = Dual.get(c).var, Series.get(order, c).var
w, x = Dual.get(b).var, Series.get(order, b).var
y, z = Dual.get(a).var, Series.get(order, a).var

def test_unary_plus():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = + y
    series = ~(+ z)
    assert abs(dual.val - y.val) < ε
    assert abs(series.val - z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_unary_minus():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = - y
    series = ~(- z)
    assert abs(dual.val + y.val) < ε
    assert abs(series.val + z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_add_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y + y
    series = ~(z + z)
    assert abs(dual.val - 2.0 * y.val) < ε
    assert abs(series.val - 2.0 * z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_add_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y + d
    series = ~(z + d)
    assert abs(dual.val - (y.val + d)) < ε
    assert abs(series.val - (z.val + d)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_add_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d + y
    series = ~(d + z)
    assert abs(dual.val - (y.val + d)) < ε
    assert abs(series.val - (z.val + d)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_add_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y + b
    series = ~(z + b)
    assert abs(dual.val - (y.val + b)) < ε
    assert abs(series.val - (z.val + b)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_add_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b + y
    series = ~(b + z)
    assert abs(dual.val - (y.val + b)) < ε
    assert abs(series.val - (z.val + b)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_subtract_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y - y
    series = ~(z - z)
    assert abs(dual.val - 0.0) < ε
    assert abs(series.val - 0.0) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_subtract_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y - d
    series = ~(z - d)
    assert abs(dual.val - (y.val - d)) < ε
    assert abs(series.val - (z.val - d)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_subtract_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d - y
    series = ~(d - z)
    assert abs(dual.val - (d - y.val)) < ε
    assert abs(series.val - (d - z.val)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_subtract_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y - b
    series = ~(z - b)
    assert abs(dual.val - (y.val - b)) < ε
    assert abs(series.val - (z.val - b)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_subtract_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b - y
    series = ~(b - z)
    assert abs(dual.val - (b - y.val)) < ε
    assert abs(series.val - (b - z.val)) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_multiply_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y * y
    series = ~(z * z)
    assert abs(dual.val - y.val**2) < ε
    assert abs(series.val - z.val**2) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_multiply_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y * d
    series = ~(z * d)
    assert abs(dual.val - y.val * d) < ε
    assert abs(series.val - z.val * d) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_multiply_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d * y
    series = ~(d * z)
    assert abs(dual.val - d * y.val) < ε
    assert abs(series.val - d * z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_multiply_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y * b
    series = ~(z * b)
    assert abs(dual.val - y.val * b) < ε
    assert abs(series.val - z.val * b) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_multiply_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b * y
    series = ~(b * z)
    assert abs(dual.val - b * y.val) < ε
    assert abs(series.val - b * z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_divide_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y / y
    series = ~(z / z)
    assert abs(dual.val - 1.0) < ε
    assert abs(series.val - 1.0) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_divide_object_int():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = y / d
    series = ~(z / d)
    assert abs(dual.val - y.val / d) < ε
    assert abs(series.val - z.val / d) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_divide_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d / y
    series = ~(d / z)
    assert abs(dual.val - d / y.val) < ε
    assert abs(series.val - d / z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_divide_object_float():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = y / b
    series = ~(z / b)
    assert abs(dual.val - y.val / b) < ε
    assert abs(series.val - z.val / b) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_divide_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b / y
    series = ~(b / z)
    assert abs(dual.val - b / y.val) < ε
    assert abs(series.val - b / z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_reciprocal():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 1.0 / y
    series = ~(1.0 / z)
    assert abs(dual.val - 1.0 / y.val) < ε
    assert abs(series.val - 1.0 / z.val) < ε
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_object_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = (y + 1)**(y - 1)
    series = ~((z + 1)**(z - 1))
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_int_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(d, int)
    dual = d**y
    series = ~d**z
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_float_object():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    assert isinstance(b, float)
    dual = b**y
    series = ~b**z
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_int():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(d, int)
    dual = w**d
    series = ~x**d
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_int_neg():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(d, int)
    dual = w**-d
    series = ~x**-d
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_float():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(a, float)
    dual = w**a
    series = ~x**a
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_float_neg():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    assert isinstance(a, float)
    dual = w**-a
    series = ~x**-a
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w**0
    series = ~x**0
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_pow_zero_float():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w**0.0
    series = ~x**0.0
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_exp():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y.exp
    series = z.exp
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_ln():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = y.ln
    series = z.ln
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_sin():
    dual = Dual.get(pi / a).var.sin
    series = ~Series.get(order, pi / a).var.sin
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_cos():
    dual = Dual.get(pi / a).var.cos
    series = ~Series.get(order, pi / a).var.cos
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_tan():
    dual = Dual.get(pi / b).var.tan
    series = ~Series.get(order, pi / b).var.tan
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_sinh():
    dual = Dual.get(pi / a).var.sinh
    series = ~Series.get(order, pi / a).var.sinh
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_cosh():
    dual = Dual.get(pi / a).var.cosh
    series = ~Series.get(order, pi / a).var.cosh
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

def test_tanh():
    dual = Dual.get(pi / b).var.tanh
    series = ~Series.get(order, pi / b).var.tanh
    assert abs(dual.val - series.jet[0]) < ε
    assert abs(dual.der - series.jet[1]) < ε

#  Zero identities

def test_pow_neg_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = 1.0 / w**2 - w**-2.0
    series = ~(1.0 / x**2 - x**-2.0)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_pow_pos_frac_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = (w**2)**0.5 - abs(w)
    series = ~((x**2)**0.5) - abs(x)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_exp_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.exp - e**w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_ln_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.exp.ln - w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_sinh_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 0.5 * (y.exp - (-y).exp) - y.sinh
    series = ~(0.5 * (z.exp - (-z).exp) - z.sinh)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_cosh_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = 0.5 * (y.exp + (-y).exp) - y.cosh
    series = ~(0.5 * (z.exp + (-z).exp) - z.cosh)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_tan_zero():
    assert isinstance(w, Dual)
    assert isinstance(x, Series)
    dual = w.tan - w.sin / w.cos
    series = ~(x.tan - x.sin / x.cos)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_tanh_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = u.tanh - u.sinh / u.cosh
    series = ~(v.tanh - v.sinh / v.cosh)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_sin_3x_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = (3 * y).sin - 3.0 * y.sin + 4.0 * y.sin**3
    series = ~((3 * z).sin - 3.0 * z.sin + 4.0 * z.sin**3)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_cos_3x_zero():
    assert isinstance(y, Dual)
    assert isinstance(z, Series)
    dual = (3 * y).cos + 3.0 * y.cos - 4.0 * y.cos**3
    series = ~((3 * z).cos + 3.0 * z.cos - 4.0 * z.cos**3)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_sinh_3x_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = (3 * u).sinh - 3.0 * u.sinh - 4.0 * u.sinh**3
    series = ~((3 * v).sinh - 3.0 * v.sinh - 4.0 * v.sinh**3)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

def test_cosh_3x_zero():
    assert isinstance(u, Dual)
    assert isinstance(v, Series)
    dual = (3 * u).cosh + 3.0 * u.cosh - 4.0 * u.cosh**3
    series = ~((3 * v).cosh + 3.0 * v.cosh - 4.0 * v.cosh**3)
    assert abs(dual.val) < ε
    assert abs(series.jet[0]) < ε
    assert abs(dual.der) < ε
    assert abs(series.jet[1]) < ε

