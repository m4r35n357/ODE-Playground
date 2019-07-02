#!/usr/bin/env python3

from sys import argv
from gmpy2 import get_context, acos, zero

get_context().precision = 236  # Set this BEFORE importing any AD stuff!
from ad import to_mpfr, Series, Dual
from functions import *

order = int(argv[1])
assert order > 1

pi = acos(-1)
a = to_mpfr("3.0")
b = to_mpfr("4.0")
c = to_mpfr("0.5")

y = Dual.get(a).var
print(f"y = {y}", file=stderr)
z = Series.get(order, a).var
print(f"z = {z}", file=stderr)

print(" unary +", file=stderr)
print(+ y, file=stderr)
print(+ z, file=stderr)
print(" unary -", file=stderr)
print(- y, file=stderr)
print(- z, file=stderr)

print(" +", file=stderr)
print(y + b, file=stderr)
print(z + b, file=stderr)
print(b + y, file=stderr)
print(b + z, file=stderr)
print(y + y, file=stderr)
print(z + z, file=stderr)

print(" -", file=stderr)
print(y - b, file=stderr)
print(z - b, file=stderr)
print(b - y, file=stderr)
print(b - z, file=stderr)
print(y - y, file=stderr)
print(z - z, file=stderr)

print(" *", file=stderr)
print(y * b, file=stderr)
print(z * b, file=stderr)
print(b * y, file=stderr)
print(b * z, file=stderr)
print(y * y, file=stderr)
print(z * z, file=stderr)

print(" /", file=stderr)
print(y / b, file=stderr)
print(z / b, file=stderr)
print(b / y, file=stderr)
print(b / z, file=stderr)
print(y / y, file=stderr)
print(z / z, file=stderr)

print(" **b", file=stderr)
print(y**b, file=stderr)
print(z**b, file=stderr)

print(" sqr", file=stderr)
print(Dual.get(b).var.sqr, file=stderr)
print(~ Series.get(order, b).var.sqr, file=stderr)

print(" sqrt", file=stderr)
print(Dual.get(b).var.sqrt, file=stderr)
print(~ Series.get(order, b).var.sqrt, file=stderr)

print(" exp", file=stderr)
print(Dual.get(b).var.exp, file=stderr)
print(~ Series.get(order, b).var.exp, file=stderr)

print(" ln", file=stderr)
print(Dual.get(b).var.ln, file=stderr)
print(~ Series.get(order, b).var.ln, file=stderr)

si, co = Dual.get(pi / a).var.sin_cos
sine, cosine = Series.get(order, pi / a).var.sin_cos

print(f" sin(pi / {a:.1})", file=stderr)
print(si, file=stderr)
print(~ sine, file=stderr)

print(f" cos(pi / {a:.1})", file=stderr)
print(co, file=stderr)
print(~ cosine, file=stderr)

si, co = Dual.get(pi / a).var.sinh_cosh
sine, cosine = Series.get(order, pi / a).var.sinh_cosh

print(f" sinh(pi / {a:.1})", file=stderr)
print(si, file=stderr)
print(~ sine, file=stderr)

print(f" cosh(pi / {a:.1})", file=stderr)
print(co, file=stderr)
print(~ cosine, file=stderr)

print(f" tan(pi / {b:.1})", file=stderr)
print(Dual.get(pi / b).var.tan, file=stderr)
print(~ Series.get(order, pi / b).var.tan, file=stderr)

print(f" sec(pi / {b:.1})^2", file=stderr)
print(~ Series.get(order, pi / b).var.sec2, file=stderr)

print(f" asin({c:.1})", file=stderr)
print(Dual.get(c).var.asin, file=stderr)
print(~ Series.get(order, c).var.asin, file=stderr)

print(f" acos({c:.1})", file=stderr)
print(Dual.get(c).var.acos, file=stderr)
print(~ Series.get(order, c).var.acos, file=stderr)

print(f" atan({c:.1})", file=stderr)
print(Dual.get(c).var.atan, file=stderr)
print(~ Series.get(order, c).var.atan, file=stderr)

print("", file=stderr)

def schwartzschild(r, e=to_mpfr(0.962250), p_r=zero(+1), l_z=to_mpfr(4)):  # no t or phi!
    # Example: ./series_test.py 1 -8 8 1001 0 1e-12 1e-12 | ./plotMany.py 8 10 >/dev/null
    return (e**2 / (1 - 2 / r) - p_r**2 * (1 - 2 / r) - (l_z / r).sqr) / 2

def analyze(model, x0, x1, steps):
    target = zero(+1)
    for k in range(steps):
        w_x = Dual.get(x0 + k * (x1 - x0) / steps).var
        w_f = model(w_x) - target
        print(f"{w_x.val:.6e} {w_f}")
        yield w_x.val, w_f.val, w_f.der

lower = to_mpfr(argv[2])
upper = to_mpfr(argv[3])
assert upper > lower
n_steps = int(argv[4])
assert n_steps > 0

for result in analyze(composite1, lower, upper, steps=n_steps):
    pass
