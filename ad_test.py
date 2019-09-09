#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov-report html:cov_html ad_test.py -v
#  Mutation Testing
#  rm -f .mutmut-cache; mutmut --test-time-base 10.0 --paths-to-mutate ad.py run --runner 'pytest ad_test.py'
from math import pi, exp, log, sin, cos, tan, sinh, cosh, tanh, factorial
from ad import t_jet, t_horner, t_abs, t_prod, t_quot, t_pwr, t_exp, t_ln, t_sin_cos, t_tan_sec2, Series, t_asin, t_acos, t_atan
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

s_05 = Series.get(order, f05).var
s_1 = Series.get(order, f1).var
s_2 = Series.get(order, f2).var
s_3 = Series.get(order, f3).var
s_4 = Series.get(order, f4).var

d_3 = complex(f3)

data1_s = Series([1.0] * order)
for i in range(1, order):
    data1_s.jet[i] = 0.5 * i * data1_s.jet[i - 1]
data2_s = Series([1.0] * order)
for i in range(1, order):
    data2_s.jet[i] = - i * data2_s.jet[i - 1]

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

@mark.parametrize('number', [1.0, δ])
def test_t_abs_positive(number):
    delta = [number] * order
    for k in range(order):
        assert t_abs(delta, k) * number > 0.0

@mark.parametrize('number', [zero, -zero])
def test_t_abs_zero(number):
    delta = [number] * order
    for k in range(order):
        assert t_abs(delta, k) * number == 0.0

@mark.parametrize('number', [-1.0, -δ])
def test_t_abs_negative(number):
    delta = [number] * order
    for k in range(order):
        assert t_abs(delta, k) * number < 0.0

def test_exceptions_add():
    with raises(RuntimeError) as e:
        _ = s_3 + d_3
    assert "Incompatible Type: <class 'complex'>" in str(e.value)

def test_exceptions_subtract():
    with raises(RuntimeError) as e:
        _ = s_3 - d_3
    assert "Incompatible Type: <class 'complex'>" in str(e.value)

def test_exceptions_multiply():
    with raises(RuntimeError) as e:
        _ = s_3 * d_3
    assert "Incompatible Type: <class 'complex'>" in str(e.value)

def test_exceptions_divide():
    with raises(RuntimeError) as e:
        _ = s_3 / d_3
    assert "Incompatible Type: <class 'complex'>" in str(e.value)

def test_exceptions_power():
    with raises(RuntimeError) as e:
        _ = s_3**d_3
    assert "Incompatible Type: <class 'complex'>" in str(e.value)

def test_get_default():
    series = Series.get(order)
    assert len(series.jet) == order
    for term in series.jet:
        assert isinstance(term, float)
        assert term == approx(0.0)

@mark.parametrize('number', [zero, -zero, f05, -f05, i5, -i5])
def test_get(number):
    series = Series.get(order, number)
    assert len(series.jet) == order
    for term in series.jet:
        assert isinstance(term, float)
    assert series.val == approx(number)
    for term in series.jet[1:]:
        assert isinstance(term, float)
        assert term == approx(0.0)

@mark.parametrize('number, length', [(data1_s, order)])
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
    series = ~ (+ data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert result == approx(original)

def test_unary_minus():
    series = ~ (- data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert result == approx(- original)

@mark.parametrize('series', [data1_s, data2_s])
def test_abs(series):
    for a, b in zip(abs(series).jet, series.jet):
        assert a == approx(b)
    for a, b in zip(abs(- series).jet, series.jet):
        assert a == approx(b)

def test_add_object_object():
    series = data1_s + data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert result == approx(s1 + s2)

@mark.parametrize('number', [i5, f3])
def test_add_object_number(number):
    series = data1_s + number
    assert len(series.jet) == order
    assert series.val == approx(data1_s.val + number)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

@mark.parametrize('number', [i5, f3])
def test_add_number_object(number):
    series = number + data1_s
    assert len(series.jet) == order
    assert series.val == approx(number + data1_s.val)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

def test_subtract_object_object():
    series = data1_s - data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert result == approx(s1 - s2)

@mark.parametrize('number', [i5, f3])
def test_subtract_object_number(number):
    series = data1_s - number
    assert len(series.jet) == order
    assert series.val == approx(data1_s.val - number)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(original)

@mark.parametrize('number', [i5, f3])
def test_subtract_number_object(number):
    series = number - data1_s
    assert len(series.jet) == order
    assert series.val == approx(number - data1_s.val)
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert result == approx(- original)

def test_multiply_object_object():
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
    series = ~ (s_3 * number)
    assert series.val == approx(f3 * number)
    assert series.jet[1] == approx(number)
    for term in series.jet[2:]:
        assert term == approx(0.0)

@mark.parametrize('number', [i5, f3])
def test_multiply_number_object(number):
    series = ~ (number * s_3)
    assert series.val == approx(number * f3)
    assert series.jet[1] == approx(number)
    for term in series.jet[2:]:
        assert term== approx(0.0)

@mark.domain
@mark.parametrize('number', [i5, δ, - δ, - i5])
def test_divide_domain_object_good(number):
    _ = data1_s / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, -zero])
def test_divide_domain_object_bad(number):
    with raises(AssertionError):
        _ = data1_s / Series.get(order, number).var

def test_divide_object_object():
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
    _ = s_4 / number

@mark.domain
@mark.parametrize('number', [zero, - zero, 0])
def test_divide_domain_object_number_bad(number):
    with raises(AssertionError):
        _ = s_4 / number

@mark.parametrize('number', [i5, f4])
def test_divide_object_number(number):
    series = ~ (s_3 / number)
    assert series.val == approx(f3 / number)
    assert series.jet[1] == approx(1.0 / number)
    for term in series.jet[2:]:
        assert term == approx(0.0)

@mark.domain
@mark.parametrize('number', [1, δ, -δ, -1])
def test_divide_domain_number_object_good(number):
    _ = 1.0 / Series.get(order, number).var

@mark.domain
@mark.parametrize('number', [zero, - zero])
def test_divide_domain_number_object_bad(number):
    with raises(AssertionError):
        _ = 1.0 / Series.get(order, number).var

@mark.parametrize('number', [i5, f4])
def test_divide_number_object(number):
    derivative = number / f3
    series = ~ (number / s_3)
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k / f3
        assert series.jet[k] == approx(derivative)

def test_reciprocal():
    derivative = 1.0 / f3
    series = ~ (1.0 / s_3)
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k / f3
        assert series.jet[k] == approx(derivative)

@mark.parametrize('number', [1, 1.0])
def test_pow_object_neg1_number(number):
    series = ~ s_4**number
    assert series.val == approx(f4)
    assert series.jet[1] == approx(1.0)
    derivative = 1.0 / f4
    series = ~ s_4**-number
    assert series.val == approx(derivative)
    for k in range(1, order):
        derivative *= - k / s_4.val
        assert series.jet[k] == approx(derivative)

@mark.domain
@mark.parametrize('number', [1, δ, 0, zero, - zero, - δ, - 1])
def test_pow_domain_object_int_good(number):
    _ = Series.get(order, number).var**2

@mark.domain
@mark.parametrize('number', [1, δ, - δ, - 1])
def test_pow_domain_object_int_neg_good(number):
    _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [0, zero, - zero])
def test_pow_domain_object_int_neg_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2

@mark.domain
@mark.parametrize('number', [1, δ])
def test_pow_domain_object_anything_good(number):
    _ = Series.get(order, number).var**data1_s
    _ = Series.get(order, number).var**2.0
    _ = Series.get(order, number).var**-2.0

@mark.domain
@mark.parametrize('number', [0, zero, - zero, - δ, -1])
def test_pow_domain_object_anything_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var**data1_s
    with raises(AssertionError):
        _ = Series.get(order, number).var**2.0
    with raises(AssertionError):
        _ = Series.get(order, number).var**-2.0

def test_pow_object_object():
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
    series = ~ s_4**number
    assert len(series.jet) == order
    assert series.val == approx(1.0)
    for term in series.jet[1:]:
        assert term == approx(0.0)

@mark.parametrize('number', [i5, - i5, f3, - f3])
def test_pow_object_number(number):
    series = ~ s_4**number
    assert len(series.jet) == order
    assert series.val == approx(f4**number)
    assert series.jet[1] == approx(number * f4**(number - 1.0))
    series = ~ s_4**-number
    assert len(series.jet) == order
    assert series.val == approx(1.0 / f4**number)

@mark.domain
@mark.parametrize('number', [1, δ])
def test_pow_domain_number_object_good(number):
    _ = number**data1_s

@mark.domain
@mark.parametrize('number', [- δ, zero, -zero])
def test_pow_domain_number_object_bad(number):
    with raises(AssertionError):
        _ = number**data1_s

@mark.parametrize('number', [i5, f3])
def test_pow_number_object(number):
    series = ~ number**data1_s
    assert len(series.jet) == order
    assert series.val == approx(number**data1_s.val)

def test_exp():
    exponent = t_jet(order)
    t_series = s_4.exp
    for k in range(order):
        exponent[k] = t_exp(exponent, s_4.jet, k)
        assert exponent[k] == approx(t_series.jet[k])
    derivative = exp(f3)
    series = ~ s_3.exp
    assert len(series.jet) == order
    for term in series.jet:
        assert term == approx(derivative)

def test_minus_exp():
    derivative = 1.0 / exp(f3)
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
    _ = Series.get(order, number).var.ln

@mark.domain
@mark.parametrize('number', [zero, - zero, - δ, - 1])
def test_ln_domain_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.ln

def test_ln():
    logarithm = t_jet(order)
    t_series = s_3.ln
    for k in range(order):
        logarithm[k] = t_ln(logarithm, s_3.jet, k)
        assert logarithm[k] == approx(t_series.jet[k])
    derivative = 1.0 / f3
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
    series = ~ Series.get(order, pi / f3).var.sinh
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert series.jet[k] == approx(sinh(pi / f3))
        elif k % 2 == 1:
            assert series.jet[k] == approx(cosh(pi / f3))

def test_cosh():
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
    series = ~ Series.get(order, pi / f4).var.tanh
    assert len(series.jet) == order
    assert series.val == approx(tanh(pi / f4))
    assert series.jet[1] == approx((1.0 - tanh(pi / f4)**2))

@mark.domain
@mark.parametrize('number', [1.0 - δ, -1.0 + δ])
def test_asin_domain_good(number):
    _ = Series.get(order, number).var.asin

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_asin_domain_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.asin

def test_asin():
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
    _ = Series.get(order, number).var.acos

@mark.domain
@mark.parametrize('number', [3.0, 1.0, -1.0, -3.0, 3, 1, -1, -3])
def test_acos_domain_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.acos

def test_acos():
    series = s_05.acos
    acos, _ = t_jet(order), t_jet(order)
    for k in range(order):
        acos[k], _[k] = t_acos(acos, _, s_05.jet, k)
        assert acos[k] == approx(series.jet[k])
    assert series.jet[0] == approx(pi / 3.0)
    assert len(series.jet) == order

def test_atan():
    series = s_1.atan
    atan, _ = t_jet(order), t_jet(order)
    for k in range(order):
        atan[k], _[k] = t_atan(atan, _, s_1.jet, k)
        assert atan[k] == approx(series.jet[k])
    assert series.jet[0] == approx(pi / 4.0)
    assert len(series.jet) == order

def test_asinh():
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
    _ = Series.get(order, number).var.acosh

@mark.domain
@mark.parametrize('number', [zero, - zero, 1.0 - δ, 1.0])
def test_acosh_domain_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.acosh

def test_acosh():
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
    _ = Series.get(order, number).var.atanh

@mark.domain
@mark.parametrize('number', [1.0, -1.0, 1, -1])
def test_atanh_domain_bad(number):
    with raises(AssertionError):
        _ = Series.get(order, number).var.atanh

def test_atanh():
    x = s_05
    series = x.atanh
    atanh, _ = t_jet(order), t_jet(order)
    for k in range(order):
        atanh[k], _[k] = t_atan(atanh, _, s_05.jet, k, hyp=True)
        assert atanh[k] == approx(series.jet[k])
    for term in (series - 0.5 * ((1 + x) / (1 - x)).ln).jet:
        assert term == approx(0.0)

def test_var():
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
def test_diff_squares():
    for term in (data1_s**2 - data2_s**2 - (data1_s - data2_s) * (data1_s + data2_s)).jet:
        assert term == approx(0.0)
    for term in (data1_s**2.0 - data2_s**2.0 - (data1_s - data2_s) * (data1_s + data2_s)).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pow_neg_zero(series):
    for term in (1.0 / series**2 - series**-2.0).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pow_frac_zero(series):
    for term in ((series**2)**0.5 - abs(series)).jet:
        assert term == approx(0.0)
    for term in ((series**2)**-0.5 - 1.0 / abs(series)).jet:
        assert term == approx(0.0)

@mark.toplevel
def test_exp_zero():
    for term in (data1_s.exp.ln - data1_s).jet:
        assert term == approx(0.0)
    for term in ((data1_s + data2_s).exp - data1_s.exp * data2_s.exp).jet:
        assert term == approx(0.0)

@mark.toplevel
def test_ln_zero():
    for term in (data1_s.ln.exp - data1_s).jet:
        assert term == approx(0.0)
    for term in ((data1_s * data2_s).ln - (data1_s.ln + data2_s.ln)).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_zero(series):
    for term in (0.5 * (series.exp - (- series).exp) - series.sinh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_zero(series):
    for term in (0.5 * (series.exp + (- series).exp) - series.cosh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tan_zero(series):
    for term in (series.tan - series.sin / series.cos).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_zero(series):
    for term in (series.tanh - series.sinh / series.cosh).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_trig(series):
    for term in (series.cos**2 + series.sin**2 - 1.0).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_pythagoras_hyp(series):
    for term in (series.cosh**2 - series.sinh**2 - 1.0).jet:
        assert term == approx(0.0)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sin_3x_zero(series):
    for a, b in zip((3 * series).sin.jet, (3.0 * series.sin - 4.0 * series.sin**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cos_3x_zero(series):
    for a, b in zip((3 * series).cos.jet, (- 3.0 * series.cos + 4.0 * series.cos**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_3x_zero(series):
    for a, b in zip((3 * series).sinh.jet, (3.0 * series.sinh + 4.0 * series.sinh**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_3x_zero(series):
    for a, b in zip((3 * series).cosh.jet, (- 3.0 * series.cosh + 4.0 * series.cosh**3).jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sin_asin_zero(series):
    for a, b in zip(series.sin.asin.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cos_acos_zero(series):
    for a, b in zip(series.cos.acos.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tan_atan_zero(series):
    for a, b in zip(series.tan.atan.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_sinh_asinh_zero(series):
    for a, b in zip(series.sinh.asinh.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_cosh_acosh_zero(series):
    for a, b in zip(series.cosh.acosh.jet, series.jet):
        assert a == approx(b)

@mark.toplevel
@mark.parametrize('series', [data1_s, data2_s])
def test_tanh_atanh_zero(series):
    for a, b in zip(series.tanh.atanh.jet, series.jet):
        assert a == approx(b)
