#!/usr/bin/env python3

from gmpy2 import get_context, acos
from sys import stderr, argv

get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from dual import Dual, make_mpfr
from series import Series
from taylor import t_jet, to_mpfr

pi = acos(-1)
a = to_mpfr("3.0")
b = to_mpfr("4.0")
c = to_mpfr("0.5")

z = Series(t_jet(7, a), variable=True)
print("z = {}".format(a), file=stderr)

print("+, -", file=stderr)
print(+ z, file=stderr)
print(- z, file=stderr)

print("+", file=stderr)
print(z + b, file=stderr)
print(b + z, file=stderr)
print(z + z, file=stderr)

print("-", file=stderr)
print(z - b, file=stderr)
print(b - z, file=stderr)
print(z - z, file=stderr)

print("*", file=stderr)
print(z * b, file=stderr)
print(b * z, file=stderr)
print(z * z, file=stderr)

print("/", file=stderr)
print(z / b, file=stderr)
print(b / z, file=stderr)
print(z / z, file=stderr)

print("**b", file=stderr)
print(z ** b, file=stderr)

print("sqr", file=stderr)
print(~ Series(t_jet(7, b), variable=True).sqr, file=stderr)
print("sqrt", file=stderr)
print(~ Series(t_jet(7, b), variable=True).sqrt, file=stderr)

sine, cosine = Series(t_jet(7, pi / a), variable=True).sin_cos
print("sin(pi / {})".format(a), file=stderr)
print(~ sine, file=stderr)
print("cos(pi / {})".format(a), file=stderr)
print(~ cosine, file=stderr)
print("tan(pi / {})".format(b), file=stderr)
print(~ Series(t_jet(7, pi / 4.0), variable=True).tan, file=stderr)

print("asin(0.5)".format(c), file=stderr)
print(Series(t_jet(7, c), variable=True).asin, file=stderr)
print("acos(0.5)".format(c), file=stderr)
print(Series(t_jet(7, c), variable=True).acos, file=stderr)
print("atan(0.5)".format(c), file=stderr)
print(Series(t_jet(7, c), variable=True).atan, file=stderr)

print("", file=stderr)

print("z = {}".format(a), file=stderr)
z = Dual.get(a, variable=True)

print("+, -", file=stderr)
print(+ z, file=stderr)
print(- z, file=stderr)

print("+", file=stderr)
print(z + b, file=stderr)
print(b + z, file=stderr)
print(z + z, file=stderr)

print("-", file=stderr)
print(z - b, file=stderr)
print(b - z, file=stderr)
print(z - z, file=stderr)

print("*", file=stderr)
print(z * b, file=stderr)
print(b * z, file=stderr)
print(z * z, file=stderr)

print("/", file=stderr)
print(z / b, file=stderr)
print(b / z, file=stderr)
print(z / z, file=stderr)

print("**b", file=stderr)
print(z ** b, file=stderr)

print("sqr", file=stderr)
print(Dual.get(b, variable=True).sqr, file=stderr)

print("sqrt", file=stderr)
print(Dual.get(b, variable=True).sqrt, file=stderr)

print("sin(pi / {:.1f})".format(a), file=stderr)
print(Dual.get(pi / a, variable=True).sin, file=stderr)
print("cos(pi / {:.1f})".format(a), file=stderr)
print(Dual.get(pi / a, variable=True).cos, file=stderr)
print("tan(pi / {:.1f})".format(b), file=stderr)
print(Dual.get(pi / 4.0, variable=True).tan, file=stderr)

print("asin(0.5)".format(c), file=stderr)
print(Dual.get(c, variable=True).asin, file=stderr)
print("acos(0.5)".format(c), file=stderr)
print(Dual.get(c, variable=True).acos, file=stderr)
print("atan(0.5)".format(c), file=stderr)
print(Dual.get(c, variable=True).atan, file=stderr)

print("", file=stderr)


def lorentz(a, value):
    # Example: ./series_test.py 0 .001 .999 1001 0 1e-12 1e-12 | ./plotMany.py 1 10 >/dev/null
    return 1 / (1 - a.sqr).sqrt - value


def cosx_x3(a, value):
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return a.cos - a * a.sqr - value


def septic(a, value):
    # Example: ./series_test.py 2 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 50000 >/dev/null
    return (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6) - value


def composite1(a, value):
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # Example: ./models.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.exp + (a.sqr - 4.0).exp).ln - value


def composite2(a, value):
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (a.sqr + (a.exp - 4).sqr).sqrt - value


def playground(a, value):
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    # return (2 * a).sin - 2 * a.sin * a.cos - value
    # return (2 * a).cos - a.cos.sqr + a.sin.sqr - value
    # return (3 * a).sin - 3 * a.sin + 4 * a.sin * a.sin.sqr - value
    # return (3 * a).cos + 3 * a.cos - 4 * a.cos * a.cos.sqr - value
    # return a.sqr.sqrt - value
    # return a.exp.ln - value
    # return a.tan - a.sin / a.cos - value
    # return a.tanh - a.sinh / a.cosh - value
    # return a.sin.asin - value
    # return a.cos.acos - value
    # return a.tan.atan - value
    return (a + 7) / (3.0 - a)


def schwartzschild(r, e=to_mpfr(0.962250), p_r=to_mpfr(0), l_z=to_mpfr(4)):  # no t or phi!
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (e**2 / (1 - 2 / r) - p_r**2 * (1 - 2 / r) - (l_z / r).sqr) / 2


def analyze(model, n, x0, x1, steps=1000, target=0, f_tol=to_mpfr(1e-12), x_tol=to_mpfr(1e-12), max_it=100, n_max=13):
    w_x = Dual.get(x0, variable=True)
    for k in range(steps):
        w_x.val = x0 + k * (x1 - x0) / steps
        w_f = model(w_x, make_mpfr(target))
        print("{:.6e} {}".format(w_x.val, w_f))
        yield w_x.val, w_f.val, w_f.der


mode = int(argv[1])
assert mode == 0 or mode == 1 or mode == 2
x0 = make_mpfr(argv[2])
x1 = make_mpfr(argv[3])
assert x1 > x0
steps = int(argv[4])
assert steps > 0
target = make_mpfr(argv[5])
f_tol = make_mpfr(argv[6])
x_tol = make_mpfr(argv[7])

for result in analyze(composite1, mode, x0, x1, steps, target, f_tol, x_tol):
    pass
