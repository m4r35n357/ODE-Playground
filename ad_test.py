#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov-report html:cov_html ad_test.py -v
#  Mutation Testing
#  rm -f .mutmut-cache; mutmut --test-time-base 10.0 --paths-to-mutate ad.py run --runner 'pytest ad_test.py'
from collections import namedtuple
from math import exp, factorial
from pytest import mark, raises, approx
from ad import Components, t_jet, t_horner, t_const, t_abs, t_pwr, Series, Dual

order = 12
# noinspection NonAsciiCharacters
ε = 1.0e-12  # small error
# noinspection NonAsciiCharacters
δ = 1.0e-6  # small amount
zero = 0.0
f1 = 1.0
f2 = 2.0
f3 = 3.0
f4 = 4.0
f05 = 0.5
i5 = 5

d_4, s_4 = Dual(f4).var, Series.get(order, f4).var
d_3, s_3 = Dual(f3).var, Series.get(order, f3).var

D0 = Dual(0.0, 0.0)
D1 = Dual(1.0, 0.0)
data1_d = Dual(1.4, -6.6)
data2_d = Dual(0.5, 7.0)

S0 = Series([0.0] * order)
S1 = Series([1.0] + [0.0] * (order - 1))
data1_s = Series([0.5] * order)
data2_s = Series([0.9] * order)
for i in range(1, order):
    data1_s.jet[i] = data1_s.jet[i - 1] / (i * i)
    data2_s.jet[i] = - data2_s.jet[i - 1] / (i * i)

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def simple_get_p():
    return Parameters(a=1.0, b=0.0, c=-1.0)

def simple_ode(x, y, z, p, k):
    return Components(x=p.a * x[k],
                      y=p.b * y[k],
                      z=p.c * z[k])

def test_t_jet_default():
    jet = t_jet(order)
    assert len(jet) == order
    for term in jet:
        assert isinstance(term, float)
        assert term == approx(0.0)

@mark.parametrize('number', [zero, -zero, f05, -f05, i5, -i5])
def test_t_jet(number):
    jet = t_jet(order, number)
    assert len(jet) == order
    assert isinstance(jet[0], float)
    assert jet[0] == approx(number)
    for term in jet[1:]:
        assert isinstance(term, float)
        assert term == approx(0.0)

def test_horner():
    assert t_horner([1, 3, 0, 2], 2) == 23
    assert t_horner([-19, 7, -4, 6], 3) == 128
    assert t_horner([3, -1, 2, -4, 0, 1], 3) == 153
    assert t_horner([1, -4, 0, 0, 2, 3, 0, -2], -2) == 201

def test_t_const():
    value = 1.23456
    ref = Series.get(order, value)
    for k in range(order):
        assert t_const(value, k) == ref.jet[k]

@mark.parametrize('number', [1.0, δ])
def test_t_abs_positive(number):
    delta = [number] + [0.0] * (order - 1)
    assert t_abs(delta, 0) > 0.0
    assert t_abs(delta, 0) == number

@mark.parametrize('number', [zero, -zero])
def test_t_abs_zero(number):
    delta = [number] + [0.0] * (order - 1)
    assert t_abs(delta, 0) == 0.0

@mark.parametrize('number', [-1.0, -δ])
def test_t_abs_negative(number):
    delta = [number] + [0.0] * (order - 1)
    assert t_abs(delta, 0) > 0.0
    assert t_abs(delta, 0) == - number

def test_exceptions_add():
    with raises(RuntimeError) as e:
        _ = d_3 + s_3
    assert "Incompatible Type: <class 'ad.Series'>" in str(e.value)
    with raises(RuntimeError) as e:
        _ = s_3 + d_3
    assert "Incompatible Type: <class 'ad.Dual'>" in str(e.value)

def test_exceptions_subtract():
    with raises(RuntimeError) as e:
        _ = d_3 - s_3
    assert "Incompatible Type: <class 'ad.Series'>" in str(e.value)
    with raises(RuntimeError) as e:
        _ = s_3 - d_3
    assert "Incompatible Type: <class 'ad.Dual'>" in str(e.value)

def test_exceptions_multiply():
    with raises(RuntimeError) as e:
        _ = d_3 * s_3
    assert "Incompatible Type: <class 'ad.Series'>" in str(e.value)
    with raises(RuntimeError) as e:
        _ = s_3 * d_3
    assert "Incompatible Type: <class 'ad.Dual'>" in str(e.value)

def test_exceptions_divide():
    with raises(RuntimeError) as e:
        _ = d_3 / s_3
    assert "Incompatible Type: <class 'ad.Series'>" in str(e.value)
    with raises(RuntimeError) as e:
        _ = s_3 / d_3
    assert "Incompatible Type: <class 'ad.Dual'>" in str(e.value)

def test_exceptions_power():
    with raises(RuntimeError) as e:
        _ = d_3**s_3
    assert "Incompatible Type: <class 'ad.Series'>" in str(e.value)
    with raises(RuntimeError) as e:
        _ = s_3**d_3
    assert "Incompatible Type: <class 'ad.Dual'>" in str(e.value)

def test_get_default():
    series = Series.get(order)
    assert len(series.jet) == order
    for term in series.jet:
        assert isinstance(term, float)
        assert term == approx(0.0)

@mark.parametrize('number', [zero, -zero, f05, -f05, i5, -i5])
def test_get(number):
    dual = Dual(number)
    series = Series.get(order, number)
    assert len(series.jet) == order
    for term in series.jet:
        assert isinstance(term, float)
    assert dual.val == approx(number)
    assert dual.dot == approx(0.0)
    assert series.val == approx(number)
    for term in series.jet[1:]:
        assert isinstance(term, float)
        assert term == approx(0.0)

@mark.parametrize('number, length', [(data1_d, 2), (data1_s, order)])
def test_to_str(number, length):
    entries = str.split(str(number))
    assert len(entries) == length
    for entry in entries:
        assert len(entry) == 10

def test_to_derivatives():
    coefficients = data1_s
    derivatives = ~ data1_s
    for k in range(order):
        assert derivatives.jet[k] == approx(factorial(k) * coefficients.jet[k])

def test_unary_plus():
    dual = + data1_d
    assert dual.val == approx(data1_d.val)
    assert dual.dot == approx(data1_d.dot)
    series = ~ (+ data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert result == approx(original)

def test_unary_minus():
    dual = - data1_d
    assert dual.val == approx(- data1_d.val)
    assert dual.dot == approx(- data1_d.dot)
    series = ~ (- data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert result == approx(- original)

def compare_series(result, target):
    for k in range(order):
        assert result.jet[k] == approx(target.jet[k])

def compare_dual(result, target):
    assert result.val == approx(target.val)
    assert result.dot == approx(target.dot)

@mark.parametrize('series', [data1_s, data2_s])
def test_abs_series(series):
    compare_series(abs(series), series.sqr.sqrt)

@mark.parametrize('dual', [data1_d, data2_d])
def test_abs_dual(dual):
    compare_dual(abs(dual), dual.sqr.sqrt)

def test_add_object_object():
    dual = data1_d + data2_d
    assert dual.val == approx(data1_d.val + data2_d.val)
    assert dual.dot == approx(data1_d.dot + data2_d.dot)
    series = data1_s + data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert result == approx(s1 + s2)

@mark.parametrize('number', [i5, f3])
def test_add_object_number(number):
    dual = data1_d + number
    assert dual.val == approx(data1_d.val + number)
    assert dual.dot == approx(data1_d.dot)
    series = data1_s + number
    assert len(series.jet) == order
    assert series.val == approx(data1_s.val + number)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

@mark.parametrize('number', [i5, f3])
def test_add_number_object(number):
    dual = number + data1_d
    assert dual.val == approx(number + data1_d.val)
    assert dual.dot == approx(data1_d.dot)
    series = number + data1_s
    assert len(series.jet) == order
    assert series.val == approx(number + data1_s.val)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

def test_subtract_object_object():
    dual = data1_d - data2_d
    assert dual.val == approx(data1_d.val - data2_d.val)
    assert dual.dot == approx(data1_d.dot - data2_d.dot)
    series = data1_s - data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert result == approx(s1 - s2)

@mark.parametrize('number', [i5, f3])
def test_subtract_object_number(number):
    dual = data1_d - number
    assert dual.val == approx(data1_d.val - number)
    assert dual.dot == approx(data1_d.dot)
    series = data1_s - number
    assert len(series.jet) == order
    assert series.val == approx(data1_s.val - number)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

@mark.parametrize('number', [i5, f3])
def test_subtract_number_object(number):
    dual = number - data1_d
    assert dual.val == approx(number - data1_d.val)
    assert dual.dot == approx(- data1_d.dot)
    series = number - data1_s
    assert len(series.jet) == order
    assert series.val == approx(number - data1_s.val)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(- original)

@mark.parametrize('series', [data1_s, data2_s])
def test_multiply_object_object_series(series):
    compare_series(series * series, series**2)

@mark.parametrize('dual', [data1_d, data2_d])
def test_multiply_object_object_dual(dual):
    compare_dual(dual * dual, dual**2)

@mark.parametrize('number', [i5, f3])
def test_multiply_object_number(number):
    dual = d_3 * number
    assert dual.val == approx(f3 * number)
    assert dual.dot == approx(number)
    series = ~ (s_3 * number)
    assert series.val == approx(f3 * number)
    assert series.jet[1] == approx(number)
    for term in series.jet[2:]:
        assert term == approx(0.0)

@mark.parametrize('number', [i5, f3])
def test_multiply_number_object(number):
    dual = number * d_3
    assert dual.val == approx(number * f3)
    assert dual.dot == approx(number)
    series = ~ (number * s_3)
    assert series.val == approx(number * f3)
    assert series.jet[1] == approx(number)
    for term in series.jet[2:]:
        assert term == approx(0.0)

@mark.domain
@mark.parametrize('number', [i5, δ, - δ, - i5])
def test_divide_domain_object_good(number):
    _ = data1_d / Dual(number).var
    _ = data1_s / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, -zero])
def test_divide_domain_object_bad(number):
    with raises(AssertionError):
        _ = data1_d / Dual(number).var
    with raises(AssertionError):
        _ = data1_s / Series.get(order, number).var

@mark.parametrize('series', [data1_s, data2_s])
def test_divide_object_object_series(series):
    compare_series(series.sin / series.cos, series.tan)

@mark.parametrize('dual', [data1_d, data2_d])
def test_divide_object_object_dual(dual):
    compare_dual(dual.sin / dual.cos, dual.tan)

@mark.domain
@mark.parametrize('number', [1, δ, -δ, -1])
def test_divide_domain_object_number_good(number):
    _ = d_3 / number
    _ = s_4 / number

@mark.domain
@mark.parametrize('number', [zero, - zero, 0])
def test_divide_domain_object_number_bad(number):
    with raises(AssertionError):
        _ = d_3 / number
    with raises(AssertionError):
        _ = s_4 / number

@mark.parametrize('number', [i5, f4])
def test_divide_object_number(number):
    dual = d_3 / number
    assert dual.val == approx(f3 / number)
    assert dual.dot == approx(1.0 / number)
    series = ~ (s_3 / number)
    assert series.val == approx(f3 / number)
    assert series.jet[1] == approx(1.0 / number)
    for term in series.jet[2:]:
        assert term == approx(0.0)

@mark.domain
@mark.parametrize('number', [1, δ, -δ, -1])
def test_divide_domain_number_object_good(number):
    _ = 1.0 / Dual(number).var
    _ = 1.0 / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, - zero])
def test_divide_domain_number_object_bad(number):
    with raises(AssertionError):
        _ = 1.0 / Dual(number).var
    with raises(AssertionError):
        _ = 1.0 / Series.get(order, number).var

@mark.parametrize('number', [i5, f4])
def test_divide_number_object(number):
    derivative = number / f3
    dual = number / d_3
    assert dual.val == approx(derivative)
    assert dual.dot == approx(- derivative / f3)
    derivative = number / f3
    series = ~ (number / s_3)
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k * number / f3
        assert series.jet[k] == approx(derivative)

@mark.parametrize('series', [data1_s, data2_s])
def test_divide_number_object_series(series):
    compare_series(1.0 / series.cos.sqr, 1.0 + series.tan.sqr)
    compare_series(1.0 / series.cosh.sqr, 1.0 - series.tanh.sqr)

@mark.parametrize('dual', [data1_d, data2_d])
def test_divide_number_object_dual(dual):
    compare_dual(1.0 / dual.cos.sqr, 1.0 + dual.tan.sqr)
    compare_dual(1.0 / dual.cosh.sqr, 1.0 - dual.tanh.sqr)

@mark.parametrize('number', [1, 1.0])
def test_pow_object_neg1_number(number):
    dual = d_4**number
    assert dual.val == approx(f4)
    assert dual.dot == approx(1.0)
    series = ~ s_4**number
    assert series.val == approx(f4)
    assert series.jet[1] == approx(1.0)
    derivative = 1.0 / f4
    dual = d_4**-number
    assert dual.val == approx(derivative)
    assert dual.dot == approx(- derivative / f4)
    series = ~ s_4**-number
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k / s_4.val
        assert series.jet[k] == approx(derivative)

@mark.domain
@mark.parametrize('number', [1, δ, 0, zero, - zero, - δ, - 1])
def test_pow_domain_object_int_good(number):
    _ = Dual(number).var**2
    _ = Series.get(order, number).var**2

@mark.domain
@mark.parametrize('number', [1, δ, - δ, - 1])
def test_pow_domain_object_int_neg_good(number):
    _ = Dual(number).var**-2
    _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [0, zero, - zero])
def test_pow_domain_object_int_neg_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var**-2
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [1, δ])
def test_pow_domain_object_anything_good(number):
    _ = Dual(number).var**data1_d
    _ = Series.get(order, number).var**data1_s
    _ = Dual(number).var**2.0
    _ = Series.get(order, number).var**2.0
    _ = Dual(number).var**-2.0
    _ = Series.get(order, number).var**-2.0

@mark.domain
@mark.parametrize('number', [0, zero, - zero, - δ, -1])
def test_pow_domain_object_anything_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var**data1_d
    with raises(AssertionError):
        _ = Series.get(order, number).var**data1_s
    with raises(AssertionError):
        _ = Dual(number).var**2.0
    with raises(AssertionError):
        _ = Series.get(order, number).var**2.0
    with raises(AssertionError):
        _ = Dual(number).var**-2.0
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2.0

def test_pow_object_object():
    dual = d_3**d_4
    assert dual.val == approx(f3**f4)
    series = ~ (s_3**s_4)
    assert len(series.jet) == order
    assert series.val == approx(f3**f4)

@mark.parametrize('number', [0, zero, -zero])
def test_pow_object_zero_dual(number):
    compare_dual(data1_d**number, D1)
    compare_dual(data2_d**number, D1)

@mark.parametrize('number', [0, zero, -zero])
def test_pow_object_zero_series(number):
    compare_series(data1_s**number, S1)
    compare_series(data2_s**number, S1)

@mark.parametrize('number', [100, -100, 100.0, -100.0, 0, -0, i5, - i5, f3, - f3])
def test_pow_object_number(number):
    dual = d_4**number
    assert dual.val == approx(f4**number)
    assert dual.dot == approx(number * f4**(number - 1.0))
    series = ~ s_4**number
    assert len(series.jet) == order
    assert series.val == approx(f4**number)
    assert series.jet[1] == approx(number * f4**(number - 1.0))
    dual = d_4**-number
    assert dual.val == approx(1.0 / f4**number)
    series = ~ s_4**-number
    assert len(series.jet) == order
    assert series.val == approx(1.0 / f4**number)

@mark.domain
@mark.parametrize('number', [1, δ])
def test_pow_domain_number_object_good(number):
    _ = number**data1_d
    _ = number**data1_s

@mark.domain
@mark.parametrize('number', [- δ, zero, -zero])
def test_pow_domain_number_object_bad(number):
    with raises(AssertionError):
        _ = number**data1_d
    with raises(AssertionError):
        _ = number**data1_s

@mark.parametrize('number', [i5, f3])
def test_pow_number_object(number):
    dual = number**data1_d
    assert dual.val == approx(number**data1_d.val)
    series = ~number**data1_s
    assert len(series.jet) == order
    assert series.val == approx(number**data1_s.val)

def test_exp_dual():
    arg = 1.0
    compare_dual(Dual(arg).var.exp, Dual(exp(arg), exp(arg)))
    compare_dual((data1_d + data2_d).exp - data1_d.exp * data2_d.exp, D0)

def test_exp_series():
    arg = 1.0
    compare_series(~ Series.get(order, arg).var.exp, Series([exp(arg)] * order))
    compare_series((data1_s + data2_s).exp - data1_s.exp * data2_s.exp, S0)

@mark.domain
@mark.parametrize('number', [1, δ])
def test_ln_domain_good(number):
    _ = Dual(number).var.ln
    _ = Series.get(order, number).var.ln

@mark.domain
@mark.parametrize('number', [zero, - zero, - δ, - 1])
def test_ln_domain_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var.ln
    with raises(AssertionError):
        _ = Series.get(order, number).var.ln

@mark.parametrize('dual', [data1_d, data2_d])
def test_ln_dual(dual):
    logarithm = dual.ln
    compare_dual(dual.exp.ln, dual)
    compare_dual(dual.sqr.ln, 2 * logarithm)
    compare_dual(dual.sqr.ln, 2.0 * logarithm)
    compare_dual(dual.sqrt.ln, 0.5 * logarithm)
    compare_dual((1.0 / dual).ln, - logarithm)
    compare_dual((dual**-3).ln, - 3 * logarithm)
    compare_dual((dual**-3.0).ln, - 3.0 * logarithm)

@mark.parametrize('series', [data1_s, data2_s])
def test_ln_series(series):
    logarithm = series.ln
    compare_series(series.exp.ln, series)
    compare_series(series.sqr.ln, 2 * logarithm)
    compare_series(series.sqr.ln, 2.0 * logarithm)
    compare_series(series.sqrt.ln, 0.5 * logarithm)
    compare_series((1.0 / series).ln, - logarithm)
    compare_series((series**-3).ln, - 3 * logarithm)
    compare_series((series**-3.0).ln, - 3.0 * logarithm)

def test_ln_dual_dual():
    compare_dual((data1_d * data2_d).ln - (data1_d.ln + data2_d.ln), D0)

def test_ln_series_series():
    compare_series((data1_s * data2_s).ln - (data1_s.ln + data2_s.ln), S0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sin_cos_series(series):
    sine, cosine = series.sin_cos
    compare_series(cosine.sqr + sine.sqr, S1)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sin_cos_dual(dual):
    compare_dual(dual.cos.sqr + dual.sin.sqr, D1)

@mark.parametrize('series', [data1_s, data2_s])
def test_tan_sec2_series(series):
    tangent, secant2 = series.tan_sec2
    compare_series(secant2 - tangent.sqr, S1)

@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_cosh_series(series):
    sine, cosine = series.sinh_cosh
    compare_series(cosine.sqr - sine.sqr, S1)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_cosh_dual(dual):
    compare_dual(dual.cosh.sqr - dual.sinh.sqr, D1)

@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_series(series):
    h_tangent, h_secant2 = series.tanh_sech2
    compare_series(h_secant2 + h_tangent.sqr, S1)

@mark.parametrize('series', [data1_s, data2_s])
def test_gd_1_series(series):
    gd_1 = (abs((series.sin + 1) / series.cos)).ln
    compare_series(series.tan.asinh, gd_1)
    compare_series(series.sin.atanh, gd_1)
    compare_series(gd_1.tanh.asin, series)
    compare_series(gd_1.sinh.atan, series)

@mark.parametrize('dual', [data1_d, data2_d])
def test_gd_1_dual(dual):
    gd_1 = (abs((dual.sin + 1) / dual.cos)).ln
    compare_dual(dual.tan.asinh, gd_1)
    compare_dual(dual.sin.atanh, gd_1)
    compare_dual(gd_1.tanh.asin, dual)
    compare_dual(gd_1.sinh.atan, dual)

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_asin_domain_good(number):
    _ = Dual(number).var.asin
    _ = Series.get(order, number).var.asin

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_asin_domain_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var.asin
        _ = Series.get(order, number).var.asin

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_acos_domain_good(number):
    _ = Dual(number).var.acos
    _ = Series.get(order, number).var.acos

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_acos_domain_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var.acos
        _ = Series.get(order, number).var.acos

@mark.domain
@mark.parametrize('number', [1.0 + δ])
def test_acosh_domain_good(number):
    _ = Dual(number).var.acosh
    _ = Series.get(order, number).var.acosh

@mark.domain
@mark.parametrize('number', [zero, - zero, 1.0 - δ, 1.0])
def test_acosh_domain_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var.acosh
        _ = Series.get(order, number).var.acosh

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_atanh_domain_good(number):
    _ = Dual(number).var.atanh
    _ = Series.get(order, number).var.atanh

@mark.domain
@mark.parametrize('number', [1.0, -1.0, 1, -1])
def test_atanh_domain_bad(number):
    with raises(AssertionError):
        _ = Dual(number).var.atanh
        _ = Series.get(order, number).var.atanh

@mark.domain
@mark.parametrize('number', [1, δ])
def test_sqrt_domain_good(number):
    _ = Series.get(order, number).var.sqrt

@mark.domain
@mark.parametrize('number', [0, zero, - zero])
def test_sqrt_domain_bad_assert(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.sqrt

@mark.domain
@mark.parametrize('number', [- δ, -1])
def test_sqrt_domain_bad_value(number):
    with raises(ValueError):
        _ = Series.get(order, number).var.sqrt

def test_var():
    dual = Dual(f3).var
    assert dual.val == approx(f3)
    assert dual.dot == approx(1.0)
    series = Series.get(order, f3).var
    assert len(series.jet) == order
    assert series.val == approx(f3)
    assert series.jet[1] == approx(1.0)
    for term in series.jet[2:]:
        assert term == approx(0.0)
    length = - 1
    try:
        length = len(Series.get(1).var.jet)
        assert False
    except AssertionError:
        assert length == -1
    assert len(Series.get(2).var.jet) == 2

#  Zero identities
def test_diff_squares_dual():
    compare_dual(data1_d**2 - data2_d**2 - (data1_d - data2_d) * (data1_d + data2_d), D0)
    compare_dual(data1_d**2.0 - data2_d**2.0 - (data1_d - data2_d) * (data1_d + data2_d), D0)

def test_diff_squares_series():
    compare_series(data1_s**2 - data2_s**2 - (data1_s - data2_s) * (data1_s + data2_s), S0)
    compare_series(data1_s**2.0 - data2_s**2.0 - (data1_s - data2_s) * (data1_s + data2_s), S0)

def test_abs_identities_dual():
    compare_dual(abs(data1_d * data2_d), abs(data1_d) * abs(data2_d))
    compare_dual(abs(data1_d / data2_d), abs(data1_d) / abs(data2_d))

def test_abs_identities_series():
    compare_series(abs(data1_s * data2_s), abs(data1_s) * abs(data2_s))
    compare_series(abs(data1_s / data2_s), abs(data1_s) / abs(data2_s))

@mark.parametrize('dual', [data1_d, data2_d])
def test_pow_neg_zero_dual(dual):
    compare_dual(dual**-2.0 - 1.0 / dual**2, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_pow_neg_zero_series(series):
    compare_series(series**-2.0 - 1.0 / series**2, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_pow_frac_zero_dual(dual):
    compare_dual((dual**2)**0.5 - abs(dual), D0)
    compare_dual((dual**2)**-0.5 - 1.0 / abs(dual), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_pow_frac_zero_series(series):
    compare_series((series**2)**0.5 - abs(series), S0)
    compare_series((series**2)**-0.5 - 1.0 / abs(series), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_dual(dual):
    compare_dual(0.5 * (dual.exp - (- dual).exp) - dual.sinh, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_series(series):
    compare_series(0.5 * (series.exp - (- series).exp) - series.sinh, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_dual(dual):
    compare_dual(0.5 * (dual.exp + (- dual).exp) - dual.cosh, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_series(series):
    compare_series(0.5 * (series.exp + (- series).exp) - series.cosh, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_tan_dual(dual):
    compare_dual(dual.tan - dual.sin / dual.cos, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_tan_series(series):
    compare_series(series.tan - series.sin / series.cos, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_tanh_dual(dual):
    compare_dual(dual.tanh - dual.sinh / dual.cosh, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_series(series):
    compare_series(series.tanh - series.sinh / series.cosh, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_pythagoras_trig_dual(dual):
    compare_dual(dual.cos.sqr + dual.sin.sqr, D1)
    compare_dual(dual.cos**2 + dual.sin**2, D1)
    compare_dual(dual.cos**2.0 + dual.sin**2.0, D1)

@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_trig_series(series):
    compare_series(series.cos.sqr + series.sin.sqr, S1)
    compare_series(series.cos**2 + series.sin**2, S1)
    compare_series(series.cos**2.0 + series.sin**2.0, S1)

@mark.parametrize('dual', [data1_d, data2_d])
def test_pythagoras_hyp_dual(dual):
    compare_dual(dual.cosh.sqr - dual.sinh.sqr, D1)
    compare_dual(dual.cosh**2 - dual.sinh**2, D1)
    compare_dual(dual.cosh**2.0 - dual.sinh**2.0, D1)

@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_hyp_series(series):
    compare_series(series.cosh.sqr - series.sinh.sqr, S1)
    compare_series(series.cosh**2 - series.sinh**2, S1)
    compare_series(series.cosh**2.0 - series.sinh**2.0, S1)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sin_3x_dual(dual):
    compare_dual((3 * dual).sin - (3 * dual.sin - 4 * dual.sin**3), D0)
    compare_dual((3.0 * dual).sin - (3.0 * dual.sin - 4.0 * dual.sin**3.0), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sin_3x_series(series):
    compare_series((3 * series).sin - (3 * series.sin - 4 * series.sin**3), S0)
    compare_series((3.0 * series).sin - (3.0 * series.sin - 4.0 * series.sin**3.0), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_cos_3x_dual(dual):
    compare_dual((3 * dual).cos - (-3 * dual.cos + 4 * dual.cos**3), D0)
    compare_dual((3.0 * dual).cos - (-3.0 * dual.cos + 4.0 * dual.cos**3.0), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_cos_3x_series(series):
    compare_series((3 * series).cos - (-3 * series.cos + 4 * series.cos**3), S0)
    compare_series((3.0 * series).cos - (-3.0 * series.cos + 4.0 * series.cos**3.0), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_3x_dual(dual):
    compare_dual((3 * dual).sinh - (3 * dual.sinh + 4 * dual.sinh**3), D0)
    compare_dual((3.0 * dual).sinh - (3.0 * dual.sinh + 4.0 * dual.sinh**3.0), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_3x_series(series):
    compare_series((3 * series).sinh - (3 * series.sinh + 4 * series.sinh**3), S0)
    compare_series((3.0 * series).sinh - (3.0 * series.sinh + 4.0 * series.sinh**3.0), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_3x_dual(dual):
    compare_dual((3 * dual).cosh - (-3 * dual.cosh + 4 * dual.cosh**3), D0)
    compare_dual((3.0 * dual).cosh - (-3.0 * dual.cosh + 4.0 * dual.cosh**3.0), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_3x_series(series):
    compare_series((3 * series).cosh - (-3 * series.cosh + 4 * series.cosh**3), S0)
    compare_series((3.0 * series).cosh - (-3.0 * series.cosh + 4.0 * series.cosh**3.0), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sin_asin_dual(dual):
    compare_dual(dual.sin.asin - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sin_asin_series(series):
    compare_series(series.sin.asin - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_cos_acos_dual(dual):
    compare_dual(dual.cos.acos - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_cos_acos_series(series):
    compare_series(series.cos.acos - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_tan_atan_dual(dual):
    compare_dual(dual.tan.atan - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_tan_atan_series(series):
    compare_series(series.tan.atan - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_asinh_dual(dual):
    compare_dual(dual.sinh.asinh - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_asinh_series(series):
    compare_series(series.sinh.asinh - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_acosh_dual(dual):
    compare_dual(dual.cosh.acosh - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_acosh_series(series):
    compare_series(series.cosh.acosh - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_tanh_atanh_dual(dual):
    compare_dual(dual.tanh.atanh - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_atanh_series(series):
    compare_series(series.tanh.atanh - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sqr_sqrt_dual(dual):
    compare_dual(dual.sqr.sqrt - abs(dual), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sqr_sqrt_series(series):
    compare_series(series.sqr.sqrt - abs(series), S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sqr_dual(dual):
    compare_dual(dual.sqr - dual * dual, D0)
    compare_dual(dual.sqr - dual**2, D0)
    compare_dual(dual.sqr - dual**2.0, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sqr_series(series):
    compare_series(series.sqr - series * series, S0)
    compare_series(series.sqr - series**2, S0)
    compare_series(series.sqr - series**2.0, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_sqrt_dual(dual):
    compare_dual(dual.sqrt - dual**0.5, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_sqrt_series(series):
    compare_series(series.sqrt - series**0.5, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_exp_ln_dual(dual):
    compare_dual(dual.exp.ln - dual, D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_exp_ln_series(series):
    compare_series(series.exp.ln - series, S0)

@mark.parametrize('dual', [data1_d, data2_d])
def test_abs_reflect_dual(dual):
    compare_dual(abs(- dual) - abs(dual), D0)
    compare_dual(abs(abs(dual)) - abs(dual), D0)

@mark.parametrize('series', [data1_s, data2_s])
def test_abs_reflect_series(series):
    compare_series(abs(- series) - abs(series), S0)
    compare_series(abs(abs(series)) - abs(series), S0)
