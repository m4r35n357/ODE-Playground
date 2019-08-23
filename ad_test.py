#!/usr/bin/env python3
#  Unit Testing
#  pytest -v --cov=ad --cov-report html:cov_html ad_test.py
#  Mutation Testing
#  mut.py --runner pytest --target ad.py --unit-test ad_test -c --disable-operator AOR CRP DDL CDI SDI SDL SVD -e
#  rm -f .mutmut-cache; mutmut --test-time-base 4.0 --paths-to-mutate ad.py run --runner 'pytest ad_test.py'
#  cosmic-ray init config.toml my_session.sqlite; cosmic-ray exec my_session.sqlite
#  cr-html my_session.sqlite > my_session.html
from math import pi, e, exp, log, sin, cos, tan, sinh, cosh, tanh
from ad import t_jet, t_horner, t_prod, t_quot, t_pwr, t_exp, t_ln, t_sin_cos, t_tan_sec2, Series, Dual

order = 6
ε = 1.0e-12  # small error
δ = 1.0e-6  # small amount
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
    assert len(jet) == order
    for term in jet:
        assert isinstance(term, float)
        assert abs(term) < ε
    jet = t_jet(order, c)
    assert len(jet) == order
    assert isinstance(jet[0], float)
    assert abs(jet[0] - c) < ε
    for term in jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε
    jet = t_jet(order, d)
    assert len(jet) == order
    assert isinstance(jet[0], float)
    assert abs(jet[0] - d) < ε
    for term in jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε

def test_horner():
    assert abs(t_horner([-19, 7, -4, 6], 3) - 128) < ε
    assert abs(t_horner([-19.0, 7.0, -4.0, 6.0], 3.0) - 128.0) < ε

def test_exceptions_add():
    try:
        _ = y + z
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = z + y
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)

def test_exceptions_subtract():
    try:
        _ = y - z
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = z - y
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_multiply():
    try:
        _ = y * z
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = z * y
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_divide():
    try:
        _ = y / z
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = z / y
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)


def test_exceptions_power():
    try:
        _ = y**z
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Series'>" in str(error)
    try:
        _ = z**y
        assert False
    except RuntimeError as error:
        assert "Incompatible Type: <class 'ad.Dual'>" in str(error)

def test_get():
    dual = Dual.get()
    series = Series.get(order)
    assert len(series.jet) == order
    assert isinstance(dual.val, float)
    assert isinstance(dual.der, float)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert isinstance(term, float)
        assert abs(term) < ε
    integer = 1
    dual = Dual.get(integer)
    series = Series.get(order, integer)
    assert len(series.jet) == order
    assert isinstance(dual.val, float)
    assert isinstance(dual.der, float)
    for term in series.jet:
        assert isinstance(term, float)
    assert abs(dual.val - float(integer)) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - float(integer)) < ε
    for term in series.jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε
    real = 1.0
    dual = Dual.get(real)
    series = Series.get(order, real)
    assert len(series.jet) == order
    assert isinstance(dual.val, float)
    assert isinstance(dual.der, float)
    for term in series.jet:
        assert isinstance(term, float)
    assert abs(dual.val - real) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - real) < ε
    for term in series.jet[1:]:
        assert isinstance(term, float)
        assert abs(term) < ε

def test_to_str():
    entries = str.split(str(y))
    assert len(entries) == 2
    for entry in entries:
        assert len(entry) == 10
    entries = str.split(str(z))
    assert len(entries) == order
    for entry in entries:
        assert len(entry) == 10

def test_unary_plus():
    dual = + y
    series = ~(+ z)
    assert len(series.jet) == order
    assert abs(dual.val - a) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - a) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_unary_minus():
    dual = - y
    series = ~(- z)
    assert len(series.jet) == order
    assert abs(dual.val + a) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val + a) < ε
    assert abs(series.jet[1] + 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_abs():
    dual = abs(u)
    series = ~(abs(v))
    assert len(series.jet) == order
    assert abs(dual.val - c) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - c) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    dual = abs(Dual.get(- c).var)
    series = ~(abs(Series.get(order, - c).var))
    assert len(series.jet) == order
    assert abs(dual.val - c) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val - c) < ε
    assert abs(series.jet[1] + 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε
    dual = abs(zero_d)
    series = ~(abs(zero_s))
    assert len(series.jet) == order
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_add_object_object():
    dual = y + w
    series = ~(z + x)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 2.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 2.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_add_object_int():
    dual = y + d
    series = ~(z + d)
    assert abs(dual.val - (a + d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + d)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_add_int_object():
    dual = d + y
    series = ~(d + z)
    assert abs(dual.val - (a + d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + d)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_add_object_float():
    dual = y + b
    series = ~(z + b)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_add_float_object():
    dual = b + y
    series = ~(b + z)
    assert abs(dual.val - (a + b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a + b)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_subtract_object_object():
    dual = y - w
    series = ~(z - x)
    assert abs(dual.val - (a - b)) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - (a - b)) < ε
    assert abs(series.jet[1]) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_subtract_object_int():
    dual = y - d
    series = ~(z - d)
    assert abs(dual.val - (a - d)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a - d)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_subtract_int_object():
    dual = d - y
    series = ~(d - z)
    assert abs(dual.val - (d - a)) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val - (d - a)) < ε
    assert abs(series.jet[1] + 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_subtract_object_float():
    dual = y - b
    series = ~(z - b)
    assert abs(dual.val - (a - b)) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - (a - b)) < ε
    assert abs(series.jet[1] - 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_subtract_float_object():
    dual = b - y
    series = ~(b - z)
    assert abs(dual.val - (b - a)) < ε
    assert abs(dual.der + 1.0) < ε
    assert abs(series.val - (b - a)) < ε
    assert abs(series.jet[1] + 1.0) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_multiply_object_object():
    t_series = z * x
    for k in range(order):
        assert t_prod(z.jet, x.jet, k) == t_series.jet[k]
    dual = y * w
    series = ~(z * x)
    assert abs(dual.val - a * b) < ε
    assert abs(dual.der - (a + b)) < ε
    assert abs(series.val - a * b) < ε
    assert abs(series.jet[1] - (a + b)) < ε
    assert abs(series.jet[2] - 2.0) < ε
    for term in series.jet[3:]:
        assert abs(term) < ε

def test_multiply_object_int():
    dual = y * d
    series = ~(z * d)
    assert abs(dual.val - a * d) < ε
    assert abs(dual.der - d) < ε
    assert abs(series.val - a * d) < ε
    assert abs(series.jet[1] - d) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_multiply_int_object():
    dual = d * y
    series = ~(d * z)
    assert abs(dual.val - d * a) < ε
    assert abs(dual.der - d) < ε
    assert abs(series.val - d * a) < ε
    assert abs(series.jet[1] - d) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_multiply_object_float():
    dual = y * b
    series = ~(z * b)
    assert abs(dual.val - a * b) < ε
    assert abs(dual.der - b) < ε
    assert abs(series.val - a * b) < ε
    assert abs(series.jet[1] - b) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_multiply_float_object():
    dual = b * y
    series = ~(b * z)
    assert abs(dual.val - b * a) < ε
    assert abs(dual.der - b) < ε
    assert abs(series.val - b * a) < ε
    assert abs(series.jet[1] - b) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_divide_domain_object_object():
    try:
        _ = y / Dual.get(δ).var
    except AssertionError:
        assert False
    try:
        _ = z / Series.get(order, δ).var
    except AssertionError:
        assert False
    try:
        _ = y / Dual.get().var
        assert False
    except AssertionError:
        pass
    try:
        _ = z / Series.get(order).var
        assert False
    except AssertionError:
        pass
    try:
        _ = y / Dual.get(- δ).var
    except AssertionError:
        assert False
    try:
        _ = z / Series.get(order, - δ).var
    except AssertionError:
        assert False

def test_divide_object_object():
    quotient = t_jet(order)
    t_series = z / x
    for k in range(order):
        quotient[k] = t_quot(quotient, z.jet, x.jet, k)
        assert quotient[k] == t_series.jet[k]
    dual = y / w
    series = ~(z / x)
    assert abs(dual.val - a / b) < ε
    assert abs(dual.der - (b - a) / b**2) < ε
    assert abs(series.val - a / b) < ε
    assert abs(series.jet[1] - (b - a) / b**2) < ε

def test_divide_object_int():
    dual = y / d
    series = ~(z / d)
    assert abs(dual.val - a / d) < ε
    assert abs(dual.der - 1.0 / d) < ε
    assert abs(series.val - a / d) < ε
    assert abs(series.jet[1] - 1.0 / d) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_divide_int_object():
    derivative = d / a
    dual = d / y
    series = ~(d / z)
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / a) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, series.n):
        derivative *= - k / a
        assert abs(series.jet[k] - derivative) < ε

def test_divide_domain_object_float():
    try:
        _ = x / (- δ)
    except AssertionError:
        assert False
    try:
        _ = y / (- δ)
    except AssertionError:
        assert False
    try:
        _ = x / 0.0
        assert False
    except AssertionError:
        pass
    try:
        _ = y / 0.0
        assert False
    except AssertionError:
        pass
    try:
        _ = x / δ
    except AssertionError:
        assert False
    try:
        _ = y / δ
    except AssertionError:
        assert False

def test_divide_object_float():
    dual = y / b
    series = ~(z / b)
    assert abs(dual.val - a / b) < ε
    assert abs(dual.der - 1.0 / b) < ε
    assert abs(series.val - a / b) < ε
    assert abs(series.jet[1] - 1.0 / b) < ε
    for term in series.jet[2:]:
        assert abs(term) < ε

def test_divide_domain_float_object():
    try:
        _ = 1.0 / Dual.get(- δ).var
    except AssertionError:
        assert False
    try:
        _ = 1.0 / Series.get(order, - δ).var
    except AssertionError:
        assert False
    try:
        _ = 1.0 / Dual.get().var
        assert False
    except AssertionError:
        pass
    try:
        _ = 1.0 / Series.get(order).var
        assert False
    except AssertionError:
        pass
    try:
        _ = 1.0 / Dual.get(δ).var
    except AssertionError:
        assert False
    try:
        _ = 1.0 / Series.get(order, δ).var
    except AssertionError:
        assert False

def test_divide_float_object():
    derivative = b / a
    dual = b / y
    series = ~(b / z)
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / a) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, series.n):
        derivative *= - k / a
        assert abs(series.jet[k] - derivative) < ε

def test_reciprocal():
    derivative = 1.0 / a
    dual = 1.0 / y
    series = ~(1.0 / z)
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / a) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, series.n):
        derivative *= - k / a
        assert abs(series.jet[k] - derivative) < ε

def test_pow_object_neg1_int():
    dual = w**1
    series = ~x**1
    assert abs(dual.val - b) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - b) < ε
    assert abs(series.jet[1] - 1.0) < ε
    derivative = 1.0 / b
    dual = w**-1
    series = ~x**-1
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / b) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, series.n):
        derivative *= - k / b
        assert abs(series.jet[k] - derivative) < ε

def test_pow_object_neg1_float():
    dual = w**1.0
    series = ~x**1.0
    assert abs(dual.val - b) < ε
    assert abs(dual.der - 1.0) < ε
    assert abs(series.val - b) < ε
    assert abs(series.jet[1] - 1.0) < ε
    derivative = 1.0 / b
    dual = w**-1.0
    series = ~x**-1.0
    assert abs(dual.val - derivative) < ε
    assert abs(dual.der + derivative / b) < ε
    assert abs(series.val - derivative) < ε
    for k in range(1, series.n):
        derivative *= - k / b
        assert abs(series.jet[k] - derivative) < ε

def test_pow_domain_object_object():
    try:
        _ = Dual.get().var**y
        assert False
    except AssertionError:
        pass
    try:
        _ = Series.get(order).var**z
        assert False
    except AssertionError:
        pass
    try:
        _ = y**Dual.get(- δ).var
        assert False
    except AssertionError:
        pass
    try:
        _ = z**Series.get(order, - δ).var
        assert False
    except AssertionError:
        pass
    try:
        _ = y**Dual.get().var
        assert False
    except AssertionError:
        pass
    try:
        _ = z**Series.get(order).var
        assert False
    except AssertionError:
        pass
    try:
        _ = y**Dual.get(δ).var
    except AssertionError:
        assert False
    try:
        _ = z**Series.get(order, δ).var
    except AssertionError:
        assert False

def test_pow_object_object():
    dual = y**w
    series = ~(z**x)
    assert abs(dual.val - a**b) < ε
    assert abs(series.val - a**b) < ε

def test_pow_int_object():
    dual = d**y
    series = ~d**z
    assert abs(dual.val - d**a) < ε
    assert abs(series.val - d**a) < ε

def test_pow_domain_float_object():
    try:
        _ = (- δ)**y
        assert False
    except AssertionError:
        pass
    try:
        _ = (- δ)**z
        assert False
    except AssertionError:
        pass
    try:
        _ = 0.0**y
        assert False
    except AssertionError:
        pass
    try:
        _ = 0.0**z
        assert False
    except AssertionError:
        pass
    try:
        _ = δ**y
    except AssertionError:
        assert False
    try:
        _ = δ**z
    except AssertionError:
        assert False

def test_pow_float_object():
    dual = b**y
    series = ~b**z
    assert abs(dual.val - b**a) < ε
    assert abs(series.val - b**a) < ε

def test_pow_domain_object_int():
    try:
        _ = Dual.get(- δ).var**2
    except AssertionError:
        assert False
    try:
        _ = Series.get(order, - δ).var**2
    except AssertionError:
        assert False
    try:
        _ = Dual.get().var**2
    except AssertionError:
        assert False
    try:
        _ = Series.get(order).var**2
    except AssertionError:
        assert False
    try:
        _ = Dual.get(δ).var**2
    except AssertionError:
        assert False
    try:
        _ = Series.get(order, δ).var**2
    except AssertionError:
        assert False

def test_pow_object_int():
    dual = w**d
    series = ~x**d
    assert abs(dual.val - b**d) < ε
    assert abs(dual.der - d * b**(d - 1)) < ε
    assert abs(series.val - b**d) < ε
    assert abs(series.jet[1] - d * b**(d - 1)) < ε
    dual = w**0
    series = ~x**0
    assert abs(dual.val - 1.0) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - 1.0) < ε
    for term in series.jet[1:]:
        assert abs(term) < ε
    dual = w**-d
    series = ~x**-d
    assert abs(dual.val - 1.0 / b**d) < ε
    assert abs(series.val - 1.0 / b**d) < ε

def test_pow_domain_object_float():
    try:
        _ = Dual.get(- δ).var**2.0
        assert False
    except AssertionError:
        pass
    try:
        _ = Series.get(order, - δ).var**2.0
        assert False
    except AssertionError:
        pass
    try:
        _ = Dual.get().var**2.0
        assert False
    except AssertionError:
        pass
    try:
        _ = Series.get(order).var**2.0
        assert False
    except AssertionError:
        pass
    try:
        _ = Dual.get(δ).var**2.0
    except AssertionError:
        assert False
    try:
        _ = Series.get(order, δ).var**2.0
    except AssertionError:
        assert False

def test_pow_object_float():
    power = t_jet(order)
    t_series = x**a
    for k in range(order):
        power[k] = t_pwr(power, x.jet, a, k)
        assert power[k] == t_series.jet[k]
    dual = w**a
    series = ~x**a
    assert abs(dual.val - b**a) < ε
    assert abs(dual.der - a * b**(a - 1.0)) < ε
    assert abs(series.val - b**a) < ε
    assert abs(series.jet[1] - a * b**(a - 1.0)) < ε
    dual = w**0.0
    series = ~x**0.0
    assert abs(dual.val - 1.0) < ε
    assert abs(dual.der) < ε
    assert abs(series.val - 1.0) < ε
    for term in series.jet[1:]:
        assert abs(term) < ε
    dual = w**-a
    series = ~x**-a
    assert abs(dual.val - 1.0 / b**a) < ε
    assert abs(series.val - 1.0 / b**a) < ε

def test_exp():
    exponent = t_jet(order)
    t_series = x.exp
    for k in range(order):
        exponent[k] = t_exp(exponent, x.jet, k)
        assert exponent[k] == t_series.jet[k]
    dual = y.exp
    series = ~z.exp
    assert len(series.jet) == order
    assert abs(dual.val - exp(a)) < ε
    assert abs(dual.der - exp(a)) < ε
    for term in series.jet:
        assert abs(term - exp(a)) < ε

def test_minus_exp():
    dual = (- y).exp
    series = ~(- z).exp
    assert len(series.jet) == order
    assert abs(dual.val - 1.0 / exp(a)) < ε
    assert abs(dual.der + 1.0 / exp(a)) < ε
    for k in range(series.n):
        if k % 2 == 0:
            assert abs(series.jet[k] - 1.0 / exp(a)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] + 1.0 / exp(a)) < ε

def test_ln_domain():
    try:
        _ = Dual.get(-δ).var.ln
        assert False
    except AssertionError:
        pass
    try:
        _ = Series.get(order, -δ).var.ln
        assert False
    except AssertionError:
        pass
    try:
        _ = Dual.get().var.ln
        assert False
    except AssertionError:
        pass
    try:
        _ = Series.get(order).var.ln
        assert False
    except AssertionError:
        pass
    try:
        _ = Dual.get(δ).var.ln
    except AssertionError:
        assert False
    try:
        _ = Series.get(order, δ).var.ln
    except AssertionError:
        assert False

def test_ln():
    logarithm = t_jet(order)
    t_series = z.ln
    for k in range(order):
        logarithm[k] = t_ln(logarithm, z.jet, k)
        assert logarithm[k] == t_series.jet[k]
    derivative = 1.0 / a
    dual = y.ln
    series = ~z.ln
    assert len(series.jet) == order
    assert abs(dual.val - log(a)) < ε
    assert abs(dual.der - derivative) < ε
    assert abs(series.val - log(a)) < ε
    assert abs(series.jet[1] - derivative) < ε
    for k in range(2, series.n):
        derivative *= - (k - 1) / a
        assert abs(series.jet[k] - derivative) < ε

def test_sin():
    sine, cosine = t_jet(order), t_jet(order)
    t_series_sin = z.sin
    t_series_cos = z.cos
    for k in range(order):
        sine[k], cosine[k] = t_sin_cos(sine, cosine, z.jet, k)
        assert sine[k] == t_series_sin.jet[k]
        assert cosine[k] == t_series_cos.jet[k]
    dual = Dual.get(pi / a).var.sin
    series = ~Series.get(order, pi / a).var.sin
    assert len(series.jet) == order
    assert abs(dual.val - sin(pi / a)) < ε
    assert abs(dual.der - cos(pi / a)) < ε
    for k in range(series.n):
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
    assert len(series.jet) == order
    assert abs(dual.val - cos(pi / a)) < ε
    assert abs(dual.der + sin(pi / a)) < ε
    for k in range(series.n):
        if k % 4 == 0:
            assert abs(series.jet[k] - cos(pi / a)) < ε
        elif k % 4 == 1:
            assert abs(series.jet[k] + sin(pi / a)) < ε
        elif k % 4 == 2:
            assert abs(series.jet[k] + cos(pi / a)) < ε
        elif k % 4 == 3:
            assert abs(series.jet[k] - sin(pi / a)) < ε

def test_tan():
    tangent, secant2 = t_jet(order), t_jet(order)
    t_series_tan = z.tan
    t_series_sec2 = z.sec2
    for k in range(order):
        tangent[k], secant2[k] = t_tan_sec2(tangent, secant2, z.jet, k)
        assert tangent[k] == t_series_tan.jet[k]
        assert secant2[k] == t_series_sec2.jet[k]
    dual = Dual.get(pi / b).var.tan
    series = ~Series.get(order, pi / b).var.tan
    assert len(series.jet) == order
    assert abs(dual.val - tan(pi / b)) < ε
    assert abs(dual.der - (1.0 + tan(pi / b)**2)) < ε
    assert abs(series.val - tan(pi / b)) < ε
    assert abs(series.jet[1] - (1.0 + tan(pi / b)**2)) < ε

def test_sinh():
    h_sine, h_cosine = t_jet(order), t_jet(order)
    t_series_sinh = z.sinh
    t_series_cosh = z.cosh
    for k in range(order):
        h_sine[k], h_cosine[k] = t_sin_cos(h_sine, h_cosine, z.jet, k, hyp=True)
        assert h_sine[k] == t_series_sinh.jet[k]
        assert h_cosine[k] == t_series_cosh.jet[k]
    dual = Dual.get(pi / a).var.sinh
    series = ~Series.get(order, pi / a).var.sinh
    assert len(series.jet) == order
    assert abs(dual.val - sinh(pi / a)) < ε
    assert abs(dual.der - cosh(pi / a)) < ε
    for k in range(series.n):
        if k % 2 == 0:
            assert abs(series.jet[k] - sinh(pi / a)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] - cosh(pi / a)) < ε

def test_cosh():
    dual = Dual.get(pi / a).var.cosh
    series = ~Series.get(order, pi / a).var.cosh
    assert len(series.jet) == order
    assert abs(dual.val - cosh(pi / a)) < ε
    assert abs(dual.der - sinh(pi / a)) < ε
    for k in range(series.n):
        if k % 2 == 0:
            assert abs(series.jet[k] - cosh(pi / a)) < ε
        elif k % 2 == 1:
            assert abs(series.jet[k] - sinh(pi / a)) < ε

def test_tanh():
    h_tangent, h_secant2 = t_jet(order), t_jet(order)
    t_series_tanh = z.tanh
    t_series_sech2 = z.sech2
    for k in range(order):
        h_tangent[k], h_secant2[k] = t_tan_sec2(h_tangent, h_secant2, z.jet, k, hyp=True)
        assert h_tangent[k] == t_series_tanh.jet[k]
        assert h_secant2[k] == t_series_sech2.jet[k]
    dual = Dual.get(pi / b).var.tanh
    series = ~Series.get(order, pi / b).var.tanh
    assert len(series.jet) == order
    assert abs(dual.val - tanh(pi / b)) < ε
    assert abs(dual.der - (1.0 - tanh(pi / b)**2)) < ε
    assert abs(series.val - tanh(pi / b)) < ε
    assert abs(series.jet[1] - (1.0 - tanh(pi / b)**2)) < ε

def test_var():
    series = Series.get(order).var
    assert len(series.jet) == order

#  Zero identities

def test_pow_neg_zero():
    dual = 1.0 / w**2 - w**-2.0
    series = ~(1.0 / x**2 - x**-2.0)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_pow_pos_frac_zero():
    dual = (w**2)**0.5 - abs(w)
    series = ~((x**2)**0.5) - abs(x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_exp_zero():
    dual = w.exp - e**w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_ln_zero():
    dual = w.exp.ln - w
    series = ~(x.exp.ln - x)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_sinh_zero():
    dual = 0.5 * (y.exp - (-y).exp) - y.sinh
    series = ~(0.5 * (z.exp - (-z).exp) - z.sinh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_cosh_zero():
    dual = 0.5 * (y.exp + (-y).exp) - y.cosh
    series = ~(0.5 * (z.exp + (-z).exp) - z.cosh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_tan_zero():
    dual = w.tan - w.sin / w.cos
    series = ~(x.tan - x.sin / x.cos)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_tanh_zero():
    dual = u.tanh - u.sinh / u.cosh
    series = ~(v.tanh - v.sinh / v.cosh)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_sin_3x_zero():
    dual = (3 * y).sin - 3.0 * y.sin + 4.0 * y.sin**3
    series = ~((3 * z).sin - 3.0 * z.sin + 4.0 * z.sin**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_cos_3x_zero():
    dual = (3 * y).cos + 3.0 * y.cos - 4.0 * y.cos**3
    series = ~((3 * z).cos + 3.0 * z.cos - 4.0 * z.cos**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_sinh_3x_zero():
    dual = (3 * u).sinh - 3.0 * u.sinh - 4.0 * u.sinh**3
    series = ~((3 * v).sinh - 3.0 * v.sinh - 4.0 * v.sinh**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε

def test_cosh_3x_zero():
    dual = (3 * u).cosh + 3.0 * u.cosh - 4.0 * u.cosh**3
    series = ~((3 * v).cosh + 3.0 * v.cosh - 4.0 * v.cosh**3)
    assert abs(dual.val) < ε
    assert abs(dual.der) < ε
    for term in series.jet:
        assert abs(term) < ε
