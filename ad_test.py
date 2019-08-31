#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov-report html:cov_html ad_test.py -v
#  Mutation Testing
#  rm -f .mutmut-cache; mutmut --test-time-base 8.0 --paths-to-mutate ad.py run --runner 'pytest ad_test.py'
from math import pi, e, exp, log, sin, cos, tan, sinh, cosh, tanh, sqrt
from ad import t_jet, t_horner, t_prod, t_quot, t_pwr, t_exp, t_ln, t_sin_cos, t_tan_sec2, Series, Dual
import pytest

order = 6
ε = 1.0e-12  # small error
δ = 1.0e-6  # small amount
zero = 0.0
f3 = 3.0
f4 = 4.0
f05 = 0.5
i5 = 5

d_0, s_0 = Dual.get().var, Series.get(order).var
d_05, s_05 = Dual.get(f05).var, Series.get(order, f05).var
d_4, s_4 = Dual.get(f4).var, Series.get(order, f4).var
d_3, s_3 = Dual.get(f3).var, Series.get(order, f3).var

data1_d = Dual(3.1, -6.6)
data2_d = Dual(0.5, 7.0)
data1_s = Series([1.0] * order)
for i in range(1, order):
    data1_s.jet[i] = 0.5 * i * data1_s.jet[i - 1]
data2_s = Series([1.0] * order)
for i in range(1, order):
    data2_s.jet[i] = - i * data2_s.jet[i - 1]

def test_t_jet_no_value():
    jet = t_jet(order)
    assert len(jet) == order
    for term in jet:
        assert isinstance(term, float)
        assert abs(term) < ε

@pytest.mark.parametrize("number", [zero, -zero, f05, -f05, i5, -i5])
def test_t_jet(number):
    jet = t_jet(order, number)
    assert len(jet) == order
    assert isinstance(jet[0], float)
    assert abs(jet[0] - number) < ε
    for term in jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε

def test_horner():
    assert abs(t_horner([-19, 7, -4, 6], 3) - 128) < ε
    assert abs(t_horner([-19.0, 7.0, -4.0, 6.0], 3.0) - 128.0) < ε

def test_exceptions_add():
    try:
        _ = d_3 + s_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = s_3 + d_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)

def test_exceptions_subtract():
    try:
        _ = d_3 - s_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = s_3 - d_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_multiply():
    try:
        _ = d_3 * s_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = s_3 * d_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_divide():
    try:
        _ = d_3 / s_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = s_3 / d_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_power():
    try:
        _ = d_3**s_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = s_3**d_3
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)

@pytest.mark.parametrize("number", [zero, -zero, f05, -f05, i5, -i5])
def test_get(number):
    dual = Dual.get(number)
    series = Series.get(order, number)
    assert len(series.jet) == order
    assert isinstance(dual.val, float)
    assert isinstance(dual.der, float)
    for term in series.jet:
        assert isinstance(term, float)
    assert abs(dual.val - number) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - number) < ε
    for term in series.jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε

@pytest.mark.parametrize("number, length", [(data1_d, 2), (data1_s, order)])
def test_to_str(number, length):
    entries = str.split(str(number))
    assert len(entries) == length
    for entry in entries:
        assert len(entry) == 10

def test_unary_plus():
    dual = + data1_d
    assert abs(dual.val - data1_d.val) < ε
    assert abs(dual.der - data1_d.der) < ε
    series = ~(+ data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert abs(result - original) < ε

def test_unary_minus():
    dual = - data1_d
    assert abs(dual.val + data1_d.val) < ε
    assert abs(dual.der + data1_d.der) < ε
    series = ~(- data1_s)
    assert len(series.jet) == order
    for result, original in zip(series.jet, (~ data1_s).jet):
        assert abs(result + original) < ε

def test_abs():
    dual = abs(d_05)
    assert abs(dual.val - f05) < ε
    assert abs(dual.der - 1.0) < ε
    series = ~abs(s_05)
    assert len(series.jet) == order
    assert abs(series.val - f05) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    dual = abs(Dual.get(- f05).var)
    assert abs(dual.val - f05) < ε
    assert abs(dual.der + 1.0) < ε
    series = ~(abs(Series.get(order, - f05).var))
    assert len(series.jet) == order
    assert abs(series.val - f05) < ε
    assert abs(series.jet[1] + 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    dual = abs(- d_05)
    assert abs(dual.val - f05) < ε
    assert abs(dual.der - 1.0) < ε
    series = ~ abs(- s_05)
    assert len(series.jet) == order
    assert abs(series.val - f05) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    dual = abs(d_0)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    series = ~(abs(s_0))
    assert len(series.jet) == order
    for term in series.jet:
        assert abs(term) < ε

def test_add_object_object():
    dual = data1_d + data2_d
    assert abs(dual.val - (data1_d.val + data2_d.val)) < ε
    assert abs(dual.der - (data1_d.der + data2_d.der)) < ε
    series = data1_s + data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert abs(result - (s1 + s2)) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_add_object_number(number):
    dual = data1_d + number
    assert abs(dual.val - (data1_d.val + number)) < ε
    assert abs(dual.der - data1_d.der) < ε
    series = data1_s + number
    assert len(series.jet) == order
    assert abs(series.val - (data1_s.val + number)) < ε
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert abs(result - original) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_add_number_object(number):
    dual = number + data1_d
    assert abs(dual.val - (number + data1_d.val)) < ε
    assert abs(dual.der - data1_d.der) < ε
    series = number + data1_s
    assert len(series.jet) == order
    assert abs(series.val - (number + data1_s.val)) < ε
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert abs(result - original) < ε

def test_subtract_object_object():
    dual = data1_d + data2_d
    assert abs(dual.val - (data1_d.val + data2_d.val)) < ε
    assert abs(dual.der - (data1_d.der + data2_d.der)) < ε
    series = data1_s + data2_s
    assert len(series.jet) == order
    for result, s1, s2 in zip(series.jet, data1_s.jet, data2_s.jet):
        assert abs(result - (s1 + s2)) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_subtract_object_number(number):
    dual = data1_d - number
    assert abs(dual.val - (data1_d.val - number)) < ε
    assert abs(dual.der - data1_d.der) < ε
    series = data1_s - number
    assert len(series.jet) == order
    assert abs(series.val - (data1_s.val - number)) < ε
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert abs(result - original) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_subtract_number_object(number):
    dual = number - data1_d
    assert abs(dual.val - (number - data1_d.val)) < ε
    assert abs(dual.der + data1_d.der) < ε
    series = number - data1_s
    assert len(series.jet) == order
    assert abs(series.val - (number - data1_s.val)) < ε
    for result, original in zip(series.jet[1:], data1_s.jet[1:]):
        assert abs(result + original) < ε

def test_multiply_object_object():
    t_series = s_3 * s_4
    for k in range(order):
        assert t_prod(s_3.jet, s_4.jet, k) == t_series.jet[k]
    dual = d_3 * d_4
    series = ~(s_3 * s_4)
    assert abs(dual.val - f3 * f4) < ε
    assert abs(dual.der - (f3 + f4)) < ε
    assert abs(series.val - f3 * f4) < ε
    assert abs(series.jet[1] - (f3 + f4)) < ε
    assert abs(series.jet[2] - 2.0) < ε
    for term in series.jet[3:]:
        assert abs(term) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_multiply_object_number(number):
    dual = d_3 * number
    series = ~(s_3 * number)
    assert abs(dual.val - f3 * number) < ε
    assert abs(dual.der - number) < ε
    assert abs(series.val - f3 * number) < ε
    assert abs(series.jet[1] - number) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

@pytest.mark.parametrize("number", [i5, f3])
def test_multiply_number_object(number):
    dual = number * d_3
    series = ~(number * s_3)
    assert abs(dual.val - number * f3) < ε
    assert abs(dual.der - number) < ε
    assert abs(series.val - number * f3) < ε
    assert abs(series.jet[1] - number) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [i5, δ, - δ, - i5])
def test_divide_domain_object_good(number):
    result = False
    try:
        _ = data1_d / Dual.get(number).var
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = data1_s / Series.get(order, number).var
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [zero, -zero])
def test_divide_domain_object_bad(number):
    result = False
    try:
        _ = data1_d / Dual.get(number).var
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = data1_s / Series.get(order, number).var
        result = False
    except AssertionError:
        result = True
    assert result

def test_divide_object_object():
    quotient = t_jet(order)
    t_series = s_3 / s_4
    for k in range(order):
        quotient[k] = t_quot(quotient, s_3.jet, s_4.jet, k)
        assert quotient[k] == t_series.jet[k]
    dual = d_3 / d_4
    series = ~(s_3 / s_4)
    assert abs(dual.val - f3 / f4) < ε
    assert abs(dual.der - (f4 - f3) / f4**2) < ε
    assert abs(series.val - f3 / f4) < ε
    assert abs(series.jet[1] - (f4 - f3) / f4**2) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ, -δ, -1])
def test_divide_domain_object_number_good(number):
    result = False
    try:
        _ = s_4 / number
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = d_3 / number
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [zero, - zero, 0])
def test_divide_domain_object_number_bad(number):
    result = False
    try:
        _ = s_4 / number
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = d_3 / number
        result = False
    except AssertionError:
        result = True
    assert result

@pytest.mark.parametrize("number", [i5, f4])
def test_divide_object_number(number):
    dual = d_3 / number
    series = ~(s_3 / number)
    assert abs(dual.val - f3 / number) < ε
    assert abs(dual.der - 1.0 / number) < ε
    assert abs(series.val - f3 / number) < ε
    assert abs(series.jet[1] - 1.0 / number) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ, -δ, -1])
def test_divide_domain_number_object_good(number):
    result = False
    try:
        _ = 1.0 / Dual.get(number).var
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = 1.0 / Series.get(order, number).var
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [zero, - zero])
def test_divide_domain_number_object_bad(number):
    result = False
    try:
        _ = 1.0 / Dual.get(number).var
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = 1.0 / Series.get(order, number).var
        result = False
    except AssertionError:
        result = True
    assert result

@pytest.mark.parametrize("number", [i5, f4])
def test_divide_number_object(number):
    derivative = number / f3
    dual = number / d_3
    series = ~(number / s_3)
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / f3) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, order):
        derivative *= - k / f3
        assert abs(series.jet[k] - derivative) < ε

def test_reciprocal():
    derivative = 1.0 / f3
    dual = 1.0 / d_3
    series = ~(1.0 / s_3)
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / f3) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, order):
        derivative *= - k / f3
        assert abs(series.jet[k] - derivative) < ε

@pytest.mark.parametrize("number", [1, 1.0])
def test_pow_object_neg1_number(number):
    dual = d_4**number
    series = ~s_4**number
    assert abs(dual.val - d_4.val) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - s_4.val) < ε
    assert abs(series.jet[1] - 1.0) < ε
    derivative = 1.0 / d_4.val
    dual = d_4**-number
    series = ~s_4**-number
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / d_4.val) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, order):
        derivative *= - k / s_4.val
        assert abs(series.jet[k] - derivative) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ, 0, zero, - zero, - δ, - 1])
def test_pow_domain_object_int_good(number):
    try:
        _ = Dual.get(number).var**2
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var**2
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ, - δ, - 1])
def test_pow_domain_object_int_neg_good(number):
    try:
        _ = Dual.get(number).var**-2
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var**-2
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [0, zero, - zero])
def test_pow_domain_object_int_neg_bad(number):
    try:
        _ = Dual.get(number).var**-2
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Series.get(order, number).var**-2
        result = False
    except AssertionError:
        result = True
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ])
def test_pow_domain_object_anything_good(number):
    result = False
    try:
        _ = Dual.get(number).var**data1_d
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var**data1_s
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Dual.get(number).var**2.0
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var**2.0
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Dual.get(number).var**-2.0
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var**-2.0
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [0, zero, - zero, - δ, -1])
def test_pow_domain_object_anything_bad(number):
    result = False
    try:
        _ = Dual.get(number).var**data1_d
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Series.get(order, number).var**data1_s
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Dual.get(number).var**2.0
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Series.get(order, number).var**2.0
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Dual.get(number).var**-2.0
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Series.get(order, number).var**-2.0
        result = False
    except AssertionError:
        result = True
    assert result

def test_pow_object_object():
    dual = d_3**d_4
    assert abs(dual.val - f3**f4) < ε
    series = ~(s_3**s_4)
    assert len(series.jet) == order
    assert abs(series.val - f3**f4) < ε

@pytest.mark.parametrize("number", [0, zero, -zero])
def test_pow_object_zero(number):
    power = t_jet(order)
    t_series = s_4**number
    for k in range(order):
        power[k] = t_pwr(power, s_4.jet, number, k)
        assert power[k] == t_series.jet[k]
    dual = d_4**number
    series = ~s_4**number
    assert len(series.jet) == order
    assert abs(dual.val - 1.0) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - 1.0) < ε
    for term in series.jet[1:]:
        assert abs(term) < ε

@pytest.mark.parametrize("number", [i5, - i5, f3, - f3])
def test_pow_object_number(number):
    dual = d_4**number
    assert abs(dual.val - f4**number) < ε
    assert abs(dual.der - number * f4**(number - 1.0)) < ε
    series = ~s_4**number
    assert len(series.jet) == order
    assert abs(series.val - f4**number) < ε
    assert abs(series.jet[1] - number * f4**(number - 1.0)) < ε
    dual = d_4**-number
    assert abs(dual.val - 1.0 / f4**number) < ε
    series = ~s_4**-number
    assert len(series.jet) == order
    assert abs(series.val - 1.0 / f4**number) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ])
def test_pow_domain_number_object_good(number):
    result = False
    try:
        _ = number**data1_d
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = number**data1_s
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [- δ, zero, -zero])
def test_pow_domain_number_object_bad(number):
    result = False
    try:
        _ = number**data1_d
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = number**data1_s
        result = False
    except AssertionError:
        result = True
    assert result

@pytest.mark.parametrize("number", [i5, f3])
def test_pow_number_object(number):
    dual = number**data1_d
    assert abs(dual.val - number**data1_d.val) < ε
    series = ~number**data1_s
    assert len(series.jet) == order
    assert abs(series.val - number**data1_s.val) < ε

def test_exp():
    exponent = t_jet(order)
    t_series = s_4.exp
    for k in range(order):
        exponent[k] = t_exp(exponent, s_4.jet, k)
        assert exponent[k] == t_series.jet[k]
    dual = d_3.exp
    assert abs(dual.val - exp(f3)) < ε
    assert abs(dual.der - exp(f3)) < ε
    series = ~s_3.exp
    assert len(series.jet) == order
    for term in series.jet:
        assert abs(term - exp(f3)) < ε

def test_minus_exp():
    dual = (- d_3).exp
    series = ~(- s_3).exp
    assert len(series.jet) == order
    assert abs(dual.val - 1.0 / exp(f3)) < ε
    assert abs(dual.der + 1.0 / exp(f3)) < ε
    for k in range(order):
        if k % 2 == 0:
            assert abs(series.jet[k] - 1.0 / exp(f3)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] + 1.0 / exp(f3)) < ε

@pytest.mark.domain
@pytest.mark.parametrize("number", [1, δ])
def test_ln_domain_good(number):
    result = False
    try:
        _ = Dual.get(number).var.ln
        result = True
    except AssertionError:
        result = False
    assert result
    try:
        _ = Series.get(order, number).var.ln
        result = True
    except AssertionError:
        result = False
    assert result

@pytest.mark.domain
@pytest.mark.parametrize("number", [zero, - zero, - δ, - 1])
def test_ln_domain_bad(number):
    result = False
    try:
        _ = Dual.get(number).var.ln
        result = False
    except AssertionError:
        result = True
    assert result
    try:
        _ = Series.get(order, number).var.ln
        result = False
    except AssertionError:
        result = True
    assert result

def test_ln():
    logarithm = t_jet(order)
    t_series = s_3.ln
    for k in range(order):
        logarithm[k] = t_ln(logarithm, s_3.jet, k)
        assert logarithm[k] == t_series.jet[k]
    derivative = 1.0 / f3
    dual = d_3.ln
    assert abs(dual.val - log(f3)) < ε
    assert abs(dual.der - derivative) < ε
    series = ~s_3.ln
    assert len(series.jet) == order
    assert abs(series.val - log(f3)) < ε
    assert abs(series.jet[1] - derivative) < ε
    for k in range(2, order):
        derivative *= - (k - 1) / f3
        assert abs(series.jet[k] - derivative) < ε

def test_sin():
    sine, cosine = t_jet(order), t_jet(order)
    t_series_sin = s_3.sin
    t_series_cos = s_3.cos
    for k in range(order):
        sine[k], cosine[k] = t_sin_cos(sine, cosine, s_3.jet, k)
        assert sine[k] == t_series_sin.jet[k]
        assert cosine[k] == t_series_cos.jet[k]
    dual = Dual.get(pi / f3).var.sin
    assert abs(dual.val - sin(pi / f3)) < ε
    assert abs(dual.der - cos(pi / f3)) < ε
    series = ~Series.get(order, pi / f3).var.sin
    assert len(series.jet) == order
    for k in range(order):
        if k % 4 == 0:
            assert abs(series.jet[k] - sin(pi / f3)) < ε
        elif k % 4 == 1:
            assert abs(series.jet[k] - cos(pi / f3)) < ε
        elif k % 4 == 2:
            assert abs(series.jet[k] + sin(pi / f3)) < ε
        elif k % 4 == 3:
            assert abs(series.jet[k] + cos(pi / f3)) < ε

def test_cos():
    dual = Dual.get(pi / f3).var.cos
    assert abs(dual.val - cos(pi / f3)) < ε
    assert abs(dual.der + sin(pi / f3)) < ε
    series = ~Series.get(order, pi / f3).var.cos
    assert len(series.jet) == order
    for k in range(order):
        if k % 4 == 0:
            assert abs(series.jet[k] - cos(pi / f3)) < ε
        elif k % 4 == 1:
            assert abs(series.jet[k] + sin(pi / f3)) < ε
        elif k % 4 == 2:
            assert abs(series.jet[k] + cos(pi / f3)) < ε
        elif k % 4 == 3:
            assert abs(series.jet[k] - sin(pi / f3)) < ε

def test_tan():
    tangent, secant2 = t_jet(order), t_jet(order)
    t_series_tan = s_3.tan
    t_series_sec2 = s_3.sec2
    for k in range(order):
        tangent[k], secant2[k] = t_tan_sec2(tangent, secant2, s_3.jet, k)
        assert tangent[k] == t_series_tan.jet[k]
        assert secant2[k] == t_series_sec2.jet[k]
    dual = Dual.get(pi / f4).var.tan
    assert abs(dual.val - tan(pi / f4)) < ε
    assert abs(dual.der - (1.0 + tan(pi / f4)**2)) < ε
    series = ~Series.get(order, pi / f4).var.tan
    assert len(series.jet) == order
    assert abs(series.val - tan(pi / f4)) < ε
    assert abs(series.jet[1] - (1.0 + tan(pi / f4)**2)) < ε

def test_sinh():
    h_sine, h_cosine = t_jet(order), t_jet(order)
    t_series_sinh = s_3.sinh
    t_series_cosh = s_3.cosh
    for k in range(order):
        h_sine[k], h_cosine[k] = t_sin_cos(h_sine, h_cosine, s_3.jet, k, hyp=True)
        assert h_sine[k] == t_series_sinh.jet[k]
        assert h_cosine[k] == t_series_cosh.jet[k]
    dual = Dual.get(pi / f3).var.sinh
    assert abs(dual.val - sinh(pi / f3)) < ε
    assert abs(dual.der - cosh(pi / f3)) < ε
    series = ~Series.get(order, pi / f3).var.sinh
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert abs(series.jet[k] - sinh(pi / f3)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] - cosh(pi / f3)) < ε

def test_cosh():
    dual = Dual.get(pi / f3).var.cosh
    assert abs(dual.val - cosh(pi / f3)) < ε
    assert abs(dual.der - sinh(pi / f3)) < ε
    series = ~Series.get(order, pi / f3).var.cosh
    assert len(series.jet) == order
    for k in range(order):
        if k % 2 == 0:
            assert abs(series.jet[k] - cosh(pi / f3)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] - sinh(pi / f3)) < ε

def test_tanh():
    h_tangent, h_secant2 = t_jet(order), t_jet(order)
    t_series_tanh = s_3.tanh
    t_series_sech2 = s_3.sech2
    for k in range(order):
        h_tangent[k], h_secant2[k] = t_tan_sec2(h_tangent, h_secant2, s_3.jet, k, hyp=True)
        assert h_tangent[k] == t_series_tanh.jet[k]
        assert h_secant2[k] == t_series_sech2.jet[k]
    dual = Dual.get(pi / f4).var.tanh
    assert abs(dual.val - tanh(pi / f4)) < ε
    assert abs(dual.der - (1.0 - tanh(pi / f4)**2)) < ε
    series = ~Series.get(order, pi / f4).var.tanh
    assert len(series.jet) == order
    assert abs(series.val - tanh(pi / f4)) < ε
    assert abs(series.jet[1] - (1.0 - tanh(pi / f4)**2)) < ε

def test_var():
    dual = Dual.get(f3).var
    assert abs(dual.val - f3) < ε
    assert abs(dual.der - 1.0) < ε
    series = Series.get(order, f3).var
    assert len(series.jet) == order
    assert abs(series.val - f3) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    length = - 1
    try:
        length = len(Series.get(1).var.jet)
        assert False
    except AssertionError:
        assert length == -1
    assert len(Series.get(2).var.jet) == 2


#  Zero identities
@pytest.mark.toplevel
def test_diff_squares():
    dual = data1_d**2 - data2_d**2 - (data1_d - data2_d) * (data1_d + data2_d)
    series = data1_s**2 - data2_s**2 - (data1_s - data2_s) * (data1_s + data2_s)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
    dual = data1_d**2.0 - data2_d**2.0 - (data1_d - data2_d) * (data1_d + data2_d)
    series = data1_s**2.0 - data2_s**2.0 - (data1_s - data2_s) * (data1_s + data2_s)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_pow_neg_zero():
    dual = 1.0 / data1_d**2 - data1_d**-2.0
    series = 1.0 / data1_s**2 - data1_s**-2.0
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_pow_frac_zero():
    dual = (d_4**2)**0.5 - abs(d_4)
    series = (s_4**2)**0.5 - abs(s_4)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
    dual = (d_4**2)**-0.5 - 1.0 / abs(d_4)
    series = (s_4**2)**-0.5 - 1.0 / abs(s_4)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_exp_zero():
    dual = data1_d.exp.ln - data1_d
    series = data1_s.exp.ln - data1_s
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
    dual = (data1_d + data2_d).exp - data1_d.exp * data2_d.exp
    series = (data1_s + data2_s).exp - data1_s.exp * data2_s.exp
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_ln_zero():
    dual = data1_d.ln.exp - data1_d
    series = data1_s.ln.exp - data1_s
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
    dual = (data1_d * data2_d).ln - (data1_d.ln + data2_d.ln)
    series = (data1_s * data2_s).ln - (data1_s.ln + data2_s.ln)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_sinh_zero():
    dual = 0.5 * (data1_d.exp - (-data1_d).exp) - data1_d.sinh
    series = 0.5 * (data1_s.exp - (-data1_s).exp) - data1_s.sinh
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_cosh_zero():
    dual = 0.5 * (data1_d.exp + (-data1_d).exp) - data1_d.cosh
    series = 0.5 * (data1_s.exp + (-data1_s).exp) - data1_s.cosh
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_tan_zero():
    dual = data1_d.tan - data1_d.sin / data1_d.cos
    series = data1_s.tan - data1_s.sin / data1_s.cos
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_tanh_zero():
    dual = d_05.tanh - d_05.sinh / d_05.cosh
    series = s_05.tanh - s_05.sinh / s_05.cosh
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_sin_3x_zero():
    dual = (3 * data1_d).sin - 3.0 * data1_d.sin + 4.0 * data1_d.sin**3
    series = (3 * data1_s).sin - 3.0 * data1_s.sin + 4.0 * data1_s.sin**3
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_cos_3x_zero():
    dual = (3 * data1_d).cos + 3.0 * data1_d.cos - 4.0 * data1_d.cos**3
    series = (3 * data1_s).cos + 3.0 * data1_s.cos - 4.0 * data1_s.cos**3
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_sinh_3x_zero():
    dual = (3 * d_05).sinh - 3.0 * d_05.sinh - 4.0 * d_05.sinh**3
    series = (3 * s_05).sinh - 3.0 * s_05.sinh - 4.0 * s_05.sinh**3
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

@pytest.mark.toplevel
def test_cosh_3x_zero():
    dual = (3 * d_05).cosh + 3.0 * d_05.cosh - 4.0 * d_05.cosh**3
    series = (3 * s_05).cosh + 3.0 * s_05.cosh - 4.0 * s_05.cosh**3
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
