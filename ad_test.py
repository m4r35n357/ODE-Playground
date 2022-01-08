#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov-report html:cov_html ad_test.py -v
#  Mutation Testing
#  rm -f .mutmut-cache; mutmut --test-time-base 10.0 --paths-to-mutate ad.py run --runner 'pytest ad_test.py'
from collections import namedtuple
from math import pi, exp, log, sin, cos, tan, sinh, cosh, tanh, factorial
from ad import Components, tsm, t_jet, t_horner, t_const, t_abs, t_prod, t_quot, t_pwr, t_exp, t_ln, t_sin_cos, t_tan_sec2, Series, Dual, t_asin, t_acos, t_atan
from pytest import mark, raises, approx

order = 6
ε = 1.0e-12  # small error
δ = 1.0e-6  # small amount
zero = 0.0
f1 = 1.0
f2 = 2.0
f3 = 3.0
f4 = 4.0
f05 = 0.5
i5 = 5

d_0, s_0 = Dual.get().var, Series.get(order).var
d_05, s_05 = Dual.get(f05).var, Series.get(order, f05).var
d_4, s_4 = Dual.get(f4).var, Series.get(order, f4).var
d_3, s_3 = Dual.get(f3).var, Series.get(order, f3).var

s_1 = Series.get(order, f1).var
s_2 = Series.get(order, f2).var

data1_d = Dual(1.4, -6.6)
data2_d = Dual(0.5, 7.0)

data1_s = Series([0.5] * order)
for i in range(1, order):
    data1_s.jet[i] = 0.7 * i * data1_s.jet[i - 1]
data2_s = Series([0.9] * order)
for i in range(1, order):
    data2_s.jet[i] = - 0.5 * i * data2_s.jet[i - 1]

class Parameters(namedtuple('ParametersType', ['a', 'b', 'c'])):
    pass

def simple_get_p(order):
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
    assert t_horner([-19, 7, -4, 6], 3) == 128
    assert t_horner([-19.0, 7.0, -4.0, 6.0], 3.0) == approx(128.0)

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
    dual = Dual.get(number)
    series = Series.get(order, number)
    assert len(series.jet) == order
    assert isinstance(dual.val, float)
    assert isinstance(dual.dot, float)
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

def test_abs():
    dual = abs(d_05)
    assert dual.val == approx(f05)
    assert dual.dot == approx(1.0)
    dual = abs(Dual.get(- f05).var)
    assert dual.val == approx(f05)
    assert dual.dot == approx(- 1.0)
    dual = abs(- d_05)
    assert dual.val == approx(f05)
    assert dual.dot == approx(1.0)
    dual = abs(d_0)
    assert dual.val == approx(0.0)
    assert dual.dot == approx(1.0)
    positive = Series.get(order, 1.0).var.asinh
    negative = Series.get(order, -1.0).var.asinh
    abs_negative = abs(negative)
    for k in range(order):
        if k % 2 == 0:
            assert t_abs(negative.jet, k) == approx(positive.jet[k])
            assert abs_negative.jet[k] == approx(positive.jet[k])
        elif k % 2 == 1:
            assert t_abs(negative.jet, k) == approx(- positive.jet[k])
            assert abs_negative.jet[k] == approx(- positive.jet[k])

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

def test_multiply_object_object():
    dual = d_3 * d_4
    assert dual.val == approx(f3 * f4)
    assert dual.dot == approx(f3 + f4)
    t_series = s_3 * s_4
    for k in range(order):
        assert t_prod(s_3.jet, s_4.jet, k) == approx(t_series.jet[k])
    series = ~ (s_3 * s_4)
    assert series.val == approx(f3 * f4)
    assert series.jet[1] == approx(f3 + f4)
    assert series.jet[2] == approx(2.0)
    for term in series.jet[3:]:
        assert term == approx(0.0)

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
        assert term== approx(0.0)

@mark.domain
@mark.parametrize('number', [i5, δ, - δ, - i5])
def test_divide_domain_object_good(number):
    _ = data1_d / Dual.get(number).var
    _ = data1_s / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, -zero])
def test_divide_domain_object_bad(number):
    with raises(AssertionError):
        _ = data1_d / Dual.get(number).var
    with raises(AssertionError):
        _ = data1_s / Series.get(order, number).var

def test_divide_object_object():
    dual = d_3 / d_4
    assert dual.val - f3 / f4 == approx(0.0)
    assert dual.dot == approx((f4 - f3) / f4**2)
    quotient = t_jet(order)
    t_series = s_3 / s_4
    for k in range(order):
        quotient[k] = t_quot(quotient, s_3.jet, s_4.jet, k)
        assert quotient[k] == approx(t_series.jet[k])
    series = ~ (s_3 / s_4)
    assert series.val == approx(f3 / f4)
    assert series.jet[1] == approx((f4 - f3) / f4**2)

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
    _ = 1.0 / Dual.get(number).var
    _ = 1.0 / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, - zero])
def test_divide_domain_number_object_bad(number):
    with raises(AssertionError):
        _ = 1.0 / Dual.get(number).var
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

def test_reciprocal():
    derivative = 1.0 / f3
    dual = 1.0 / d_3
    assert dual.val == approx(derivative)
    assert dual.dot == approx(- derivative / f3)
    series = ~ (1.0 / s_3)
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k / f3
        assert series.jet[k] == approx(derivative)

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
    _ = Dual.get(number).var**2
    _ = Series.get(order, number).var**2

@mark.domain
@mark.parametrize('number', [1, δ, - δ, - 1])
def test_pow_domain_object_int_neg_good(number):
    _ = Dual.get(number).var**-2
    _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [0, zero, - zero])
def test_pow_domain_object_int_neg_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var**-2
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [1, δ])
def test_pow_domain_object_anything_good(number):
    _ = Dual.get(number).var**data1_d
    _ = Series.get(order, number).var**data1_s
    _ = Dual.get(number).var**2.0
    _ = Series.get(order, number).var**2.0
    _ = Dual.get(number).var**-2.0
    _ = Series.get(order, number).var**-2.0

@mark.domain
@mark.parametrize('number', [0, zero, - zero, - δ, -1])
def test_pow_domain_object_anything_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var**data1_d
    with raises(AssertionError):
        _ = Series.get(order, number).var**data1_s
    with raises(AssertionError):
        _ = Dual.get(number).var**2.0
    with raises(AssertionError):
        _ = Series.get(order, number).var**2.0
    with raises(AssertionError):
        _ = Dual.get(number).var**-2.0
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2.0

def test_pow_object_object():
    dual = d_3**d_4
    assert dual.val == approx(f3**f4)
    series = ~ (s_3**s_4)
    assert len(series.jet) == order
    assert series.val == approx(f3**f4)

@mark.parametrize('number', [0, zero, -zero])
def test_pow_object_zero(number):
    power = t_jet(order)
    t_series = s_4**number
    for k in range(order):
        power[k] = t_pwr(power, s_4.jet, number, k)
        assert power[k] == approx(t_series.jet[k])
    dual = d_4**number
    assert dual.val == approx(1.0)
    assert dual.dot == approx(0.0)
    series = ~ s_4**number
    assert len(series.jet) == order
    assert series.val == approx(1.0)
    for term in series.jet[1:]:
        assert term == approx(0.0)

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

def test_exp():
    exponent = t_jet(order)
    t_series = s_4.exp
    for k in range(order):
        exponent[k] = t_exp(exponent, s_4.jet, k)
        assert exponent[k] == approx(t_series.jet[k])
    derivative = exp(f3)
    dual = d_3.exp
    assert dual.val == approx(derivative)
    assert dual.dot == approx(derivative)
    series = ~ s_3.exp
    assert len(series.jet) == order
    for term in series.jet:
        assert term == approx(derivative)

def test_minus_exp():
    derivative = 1.0 / exp(f3)
    dual = (- d_3).exp
    assert dual.val == approx(derivative)
    assert dual.dot == approx(- derivative)
    series = ~ (- s_3).exp
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert series.jet[k] == approx(derivative)
        elif k % 2 == 1:
            assert series.jet[k] == approx(- derivative)

@mark.domain
@mark.parametrize('number', [1, δ])
def test_ln_domain_good(number):
    _ = Dual.get(number).var.ln
    _ = Series.get(order, number).var.ln

@mark.domain
@mark.parametrize('number', [zero, - zero, - δ, - 1])
def test_ln_domain_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var.ln
    with raises(AssertionError):
        _ = Series.get(order, number).var.ln

def test_ln():
    logarithm = t_jet(order)
    t_series = s_3.ln
    for k in range(order):
        logarithm[k] = t_ln(logarithm, s_3.jet, k)
        assert logarithm[k] == approx(t_series.jet[k])
    derivative = 1.0 / f3
    dual = d_3.ln
    assert dual.val == approx(log(f3))
    assert dual.dot == approx(derivative)
    series = ~ s_3.ln
    assert len(series.jet) == order
    assert series.val == approx(log(f3))
    assert series.jet[1] == approx(derivative)
    for k in range(2, order):
        derivative *= - (k - 1) / f3
        assert series.jet[k] == approx(derivative)

def test_sin():
    sine, cosine = t_jet(order), t_jet(order)
    t_series_sin = s_3.sin
    t_series_cos = s_3.cos
    for k in range(order):
        sine[k], cosine[k] = t_sin_cos(sine, cosine, s_3.jet, k)
        assert sine[k] == approx(t_series_sin.jet[k])
        assert cosine[k] == approx(t_series_cos.jet[k])
    dual = Dual.get(pi / f3).var.sin
    assert dual.val == approx(sin(pi / f3))
    assert dual.dot == approx(cos(pi / f3))
    series = ~ Series.get(order, pi / f3).var.sin
    assert len(series.jet) == order
    for k in range(order):
        if k % 4 == 0:
            assert series.jet[k] == approx(sin(pi / f3))
        elif k % 4 == 1:
            assert series.jet[k] == approx(cos(pi / f3))
        elif k % 4 == 2:
            assert series.jet[k] == approx(- sin(pi / f3))
        elif k % 4 == 3:
            assert series.jet[k] == approx(- cos(pi / f3))

def test_cos():
    dual = Dual.get(pi / f3).var.cos
    assert dual.val == approx(cos(pi / f3))
    assert dual.dot == approx(- sin(pi / f3))
    series = ~ Series.get(order, pi / f3).var.cos
    assert len(series.jet) == order
    for k in range(order):
        if k % 4 == 0:
            assert series.jet[k] == approx(cos(pi / f3))
        elif k % 4 == 1:
            assert series.jet[k] == approx(- sin(pi / f3))
        elif k % 4 == 2:
            assert series.jet[k] == approx(- cos(pi / f3))
        elif k % 4 == 3:
            assert series.jet[k] == approx(sin(pi / f3))

def test_tan():
    tangent, secant2 = t_jet(order), t_jet(order)
    t_series_tan = s_3.tan
    t_series_sec2 = s_3.sec2
    for k in range(order):
        tangent[k], secant2[k] = t_tan_sec2(tangent, secant2, s_3.jet, k)
        assert tangent[k] == approx(t_series_tan.jet[k])
        assert secant2[k] == approx(t_series_sec2.jet[k])
    dual = Dual.get(pi / f4).var.tan
    assert dual.val == approx(tan(pi / f4))
    assert dual.dot == approx((1.0 + tan(pi / f4) ** 2))
    series = ~ Series.get(order, pi / f4).var.tan
    assert len(series.jet) == order
    assert series.val == approx(tan(pi / f4))
    assert series.jet[1] == approx((1.0 + tan(pi / f4)**2))

def test_sinh():
    h_sine, h_cosine = t_jet(order), t_jet(order)
    t_series_sinh = s_3.sinh
    t_series_cosh = s_3.cosh
    for k in range(order):
        h_sine[k], h_cosine[k] = t_sin_cos(h_sine, h_cosine, s_3.jet, k, hyp=True)
        assert h_sine[k] == approx(t_series_sinh.jet[k])
        assert h_cosine[k] == approx(t_series_cosh.jet[k])
    dual = Dual.get(pi / f3).var.sinh
    assert dual.val == approx(sinh(pi / f3))
    assert dual.dot == approx(cosh(pi / f3))
    series = ~ Series.get(order, pi / f3).var.sinh
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert series.jet[k] == approx(sinh(pi / f3))
        elif k % 2 == 1:
            assert series.jet[k] == approx(cosh(pi / f3))

def test_cosh():
    dual = Dual.get(pi / f3).var.cosh
    assert dual.val == approx(cosh(pi / f3))
    assert dual.dot == approx(sinh(pi / f3))
    series = ~ Series.get(order, pi / f3).var.cosh
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert series.jet[k] == approx(cosh(pi / f3))
        elif k % 2 == 1:
            assert series.jet[k] == approx(sinh(pi / f3))

def test_tanh():
    h_tangent, h_secant2 = t_jet(order), t_jet(order)
    t_series_tanh = s_3.tanh
    t_series_sech2 = s_3.sech2
    for k in range(order):
        h_tangent[k], h_secant2[k] = t_tan_sec2(h_tangent, h_secant2, s_3.jet, k, hyp=True)
        assert h_tangent[k] == approx(t_series_tanh.jet[k])
        assert h_secant2[k] == approx(t_series_sech2.jet[k])
    dual = Dual.get(pi / f4).var.tanh
    assert dual.val == approx(tanh(pi / f4))
    assert dual.dot == approx((1.0 - tanh(pi / f4) ** 2))
    series = ~ Series.get(order, pi / f4).var.tanh
    assert len(series.jet) == order
    assert series.val == approx(tanh(pi / f4))
    assert series.jet[1] == approx((1.0 - tanh(pi / f4)**2))

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_asin_domain_good(number):
    _ = Dual.get(number).var.asin
    _ = Series.get(order, number).var.asin

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_asin_domain_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var.asin
        _ = Series.get(order, number).var.asin

def test_asin():
    dual = d_05.sin.asin
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    series = s_05.asin
    asin, _ = t_jet(order), t_jet(order)
    for k in range(order):
        asin[k], _[k] = t_asin(asin, _, s_05.jet, k)
        assert asin[k] == approx(series.jet[k])
    assert series.jet[0] == approx(pi / 6.0)
    assert len(series.jet) == order

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_acos_domain_good(number):
    _ = Dual.get(number).var.acos
    _ = Series.get(order, number).var.acos

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_acos_domain_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var.acos
        _ = Series.get(order, number).var.acos

def test_acos():
    dual = d_05.cos.acos
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    series = s_05.acos
    acos, _ = t_jet(order), t_jet(order)
    for k in range(order):
        acos[k], _[k] = t_acos(acos, _, s_05.jet, k)
        assert acos[k] == approx(series.jet[k])
    assert series.jet[0] == approx(pi / 3.0)
    assert len(series.jet) == order

def test_atan():
    dual = d_05.tan.atan
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    series = s_1.atan
    atan, _ = t_jet(order), t_jet(order)
    for k in range(order):
        atan[k], _[k] = t_atan(atan, _, s_1.jet, k)
        assert atan[k] == approx(series.jet[k])
    assert series.jet[0] == approx(pi / 4.0)
    assert len(series.jet) == order

def test_asinh():
    dual = d_05.sinh.asinh
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    x = s_05
    series = x.asinh
    asinh, _ = t_jet(order), t_jet(order)
    for k in range(order):
        asinh[k], _[k] = t_asin(asinh, _, s_05.jet, k, hyp=True)
        assert asinh[k] == approx(series.jet[k])
    for term in (series - (x + (x**2 + 1)**0.5).ln).jet:
        assert term == approx(0.0)

@mark.domain
@mark.parametrize('number', [1.0 + δ])
def test_acosh_domain_good(number):
    _ = Dual.get(number).var.acosh
    _ = Series.get(order, number).var.acosh

@mark.domain
@mark.parametrize('number', [zero, - zero, 1.0 - δ, 1.0])
def test_acosh_domain_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var.acosh
        _ = Series.get(order, number).var.acosh

def test_acosh():
    dual = d_05.cosh.acosh
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    x = s_2
    series = x.acosh
    acosh, _ = t_jet(order), t_jet(order)
    for k in range(order):
        acosh[k], _[k] = t_acos(acosh, _, s_2.jet, k, hyp=True)
        assert acosh[k] == approx(series.jet[k])
    for term in (series - (x + (x**2 - 1)**0.5).ln).jet:
        assert term == approx(0.0)

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_atanh_domain_good(number):
    _ = Dual.get(number).var.atanh
    _ = Series.get(order, number).var.atanh

@mark.domain
@mark.parametrize('number', [1.0, -1.0, 1, -1])
def test_atanh_domain_bad(number):
    with raises(AssertionError):
        _ = Dual.get(number).var.atanh
        _ = Series.get(order, number).var.atanh

def test_atanh():
    dual = d_05.tanh.atanh
    assert dual.val == approx(d_05.val)
    assert dual.dot == approx(d_05.dot)
    x = s_05
    series = x.atanh
    atanh, _ = t_jet(order), t_jet(order)
    for k in range(order):
        atanh[k], _[k] = t_atan(atanh, _, s_05.jet, k, hyp=True)
        assert atanh[k] == approx(series.jet[k])
    for term in (series - 0.5 * ((1 + x) / (1 - x)).ln).jet:
        assert term == approx(0.0)

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
    dual = Dual.get(f3).var
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
@mark.toplevel
def test_diff_squares_dual():
    dual = data1_d**2 - data2_d**2 - (data1_d - data2_d) * (data1_d + data2_d)
    assert dual.val == approx(0.0)
    assert dual.dot == approx(0.0)
    dual = data1_d**2.0 - data2_d**2.0 - (data1_d - data2_d) * (data1_d + data2_d)
    assert dual.val == approx(0.0)
    assert dual.dot == approx(0.0)

@mark.toplevel
def test_diff_squares_series():
    for term in (data1_s**2 - data2_s**2 - (data1_s - data2_s) * (data1_s + data2_s)).jet:
        assert term == approx(0.0)
    for term in (data1_s**2.0 - data2_s**2.0 - (data1_s - data2_s) * (data1_s + data2_s)).jet:
        assert term == approx(0.0)

@mark.toplevel
def test_abs_identities_dual():
    a, b = abs(data1_d * data2_d), abs(data1_d) * abs(data2_d)
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)
    a, b = abs(data1_d / data2_d), abs(data1_d) / abs(data2_d)
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
def test_abs_identities_series():
    data = data1_s * data2_s
    for a, b in zip(abs(data).jet, (abs(data1_s) * abs(data2_s)).jet):
        assert a == approx(b)
    data = data1_s / data2_s
    for a, b in zip(abs(data).jet, (abs(data1_s) / abs(data2_s)).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_pow_neg_zero_dual(dual):
    e = dual**-2.0 - 1.0 / dual**2
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pow_neg_zero_series(series):
    for term in (series**-2.0 - 1.0 / series**2).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_pow_frac_zero_dual(dual):
    e = (dual**2)**0.5 - abs(dual)
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)
    e = (dual**2)**-0.5 - 1.0 / abs(dual)
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pow_frac_zero_series(series):
    for term in ((series**2)**0.5 - abs(series)).jet:
        assert term == approx(0.0)
    for term in ((series**2)**-0.5 - 1.0 / abs(series)).jet:
        assert term == approx(0.0)

@mark.toplevel
def test_exp_zero_dual():
    dual = (data1_d + data2_d).exp - data1_d.exp * data2_d.exp
    assert dual.val == approx(0.0)
    assert dual.dot == approx(0.0)

@mark.toplevel
def test_exp_zero_series():
    for term in ((data1_s + data2_s).exp - data1_s.exp * data2_s.exp).jet:
        assert term == approx(0.0)

@mark.toplevel
def test_ln_zero_dual():
    dual = (data1_d * data2_d).ln - (data1_d.ln + data2_d.ln)
    assert dual.val == approx(0.0)
    assert dual.dot == approx(0.0)

@mark.toplevel
def test_ln_zero_series():
    for term in ((data1_s * data2_s).ln - (data1_s.ln + data2_s.ln)).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_zero_dual(dual):
    e = 0.5 * (dual.exp - (- dual).exp) - dual.sinh
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_zero_series(series):
    for term in (0.5 * (series.exp - (- series).exp) - series.sinh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_zero_dual(dual):
    e = 0.5 * (dual.exp + (- dual).exp) - dual.cosh
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_zero_series(series):
    for term in (0.5 * (series.exp + (- series).exp) - series.cosh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_zero(series):
    dual = 0.5 * (data1_d.exp + (-data1_d).exp) - data1_d.cosh
    assert dual.val == approx(0.0)
    assert dual.dot == approx(0.0)
    for term in (0.5 * (series.exp + (- series).exp) - series.cosh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_tan_zero_dual(dual):
    e = dual.tan - dual.sin / dual.cos
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tan_zero_series(series):
    for term in (series.tan - series.sin / series.cos).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_tanh_zero_dual(dual):
    e = dual.tanh - dual.sinh / dual.cosh
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_zero_series(series):
    for term in (series.tanh - series.sinh / series.cosh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_pythagoras_trig_dual(dual):
    e = dual.cos ** 2 + dual.sin ** 2 - 1.0
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_trig_series(series):
    for term in (series.cos**2 + series.sin**2 - 1.0).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_pythagoras_hyp_dual(dual):
    e = dual.cosh ** 2 - dual.sinh ** 2 - 1.0
    assert e.val == approx(0.0)
    assert e.dot == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_hyp_series(series):
    for term in (series.cosh**2 - series.sinh**2 - 1.0).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sin_3x_dual(dual):
    a = (3 * dual).sin
    b = 3.0 * dual.sin - 4.0 * dual.sin**3
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sin_3x_series(series):
    for a, b in zip((3 * series).sin.jet, (3.0 * series.sin - 4.0 * series.sin**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_cos_3x_dual(dual):
    a = (3 * dual).cos
    b = - 3.0 * dual.cos + 4.0 * dual.cos**3
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cos_3x_series(series):
    for a, b in zip((3 * series).cos.jet, (- 3.0 * series.cos + 4.0 * series.cos**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_3x_dual(dual):
    a = (3 * dual).sinh
    b = 3.0 * dual.sinh + 4.0 * dual.sinh**3
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_3x_series(series):
    for a, b in zip((3 * series).sinh.jet, (3.0 * series.sinh + 4.0 * series.sinh**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_3x_dual(dual):
    a = (3 * dual).cosh
    b = - 3.0 * dual.cosh + 4.0 * dual.cosh**3
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_3x_series(series):
    for a, b in zip((3 * series).cosh.jet, (- 3.0 * series.cosh + 4.0 * series.cosh**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sin_asin_dual(dual):
    a = dual.sin.asin
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sin_asin_series(series):
    for a, b in zip(series.sin.asin.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_cos_acos_dual(dual):
    a = dual.cos.acos
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cos_acos_series(series):
    for a, b in zip(series.cos.acos.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_tan_atan_dual(dual):
    a = dual.tan.atan
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tan_atan_series(series):
    for a, b in zip(series.tan.atan.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sinh_asinh_dual(dual):
    a = dual.sinh.asinh
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_asinh_series(series):
    for a, b in zip(series.sinh.asinh.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_cosh_acosh_dual(dual):
    a = dual.cosh.acosh
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_acosh_series(series):
    for a, b in zip(series.cosh.acosh.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_tanh_atanh_dual(dual):
    a = dual.tanh.atanh
    b = dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_atanh_series(series):
    for a, b in zip(series.tanh.atanh.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sqr_sqrt_dual(dual):
    a, b = dual.sqr.sqrt, abs(dual)
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sqr_sqrt_series(series):
    for a, b in zip(series.sqr.sqrt.jet, abs(series).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sqr_dual(dual):
    a, b = dual.sqr, dual**2
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sqr_series(series):
    for a, b in zip(series.sqr.jet, (series**2).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_sqrt_dual(dual):
    a, b = dual.sqrt, dual**0.5
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sqrt_series(series):
    for a, b in zip(series.sqrt.jet, (series**0.5).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('dual', [data1_d, data2_d])
def test_exp_ln_dual(dual):
    a, b = dual.exp.ln, dual
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_exp_ln_series(series):
    for a, b in zip(series.exp.ln.jet, series.jet):
        assert a == approx(b)

@mark.parametrize('dual', [data1_d, data2_d])
def test_abs_reflect_dual(dual):
    a, b = abs(- dual), abs(dual)
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)
    a, b = abs(abs(dual)), abs(dual)
    assert a.val == approx(b.val)
    assert a.dot == approx(b.dot)

@mark.parametrize('series', [data1_s, data2_s])
def test_abs_reflect_series(series):
    for a, b in zip(abs(- series).jet, abs(series).jet):
        assert a == approx(b)
    for a, b in zip(abs(abs(series)).jet, abs(series).jet):
        assert a == approx(b)
