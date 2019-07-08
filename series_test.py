#!/usr/bin/env python3

from sys import argv
from math import acos
from ad import Series, Dual
from functions import *

order = int(argv[1])
assert order > 1

pi = acos(-1)
a = 3.0
b = 4.0
c = 0.5

y = Dual.get(a).var
print(f"y = {y}", file=stderr)
z = Series.get(order, a).var
print(f"z = {z}", file=stderr)

print(" unary +", file=stderr)
print(+ y, file=stderr)
print(~(+ z), file=stderr)
print(" unary -", file=stderr)
print(- y, file=stderr)
print(~(- z), file=stderr)

print(" +", file=stderr)
print(y + b, file=stderr)
print(~(z + b), file=stderr)
print(b + y, file=stderr)
print(~(b + z), file=stderr)
print(y + y, file=stderr)
print(~(z + z), file=stderr)

print(" -", file=stderr)
print(y - b, file=stderr)
print(~(z - b), file=stderr)
print(b - y, file=stderr)
print(~(b - z), file=stderr)
print(y - y, file=stderr)
print(~(z - z), file=stderr)

print(" *", file=stderr)
print(y * b, file=stderr)
print(~(z * b), file=stderr)
print(b * y, file=stderr)
print(~(b * z), file=stderr)
print(y * y, file=stderr)
print(~(z * z), file=stderr)

print(" /", file=stderr)
print(y / b, file=stderr)
print(~(z / b), file=stderr)
print(b / y, file=stderr)
print(~(b / z), file=stderr)
print(y / y, file=stderr)
print(~(z / z), file=stderr)

print(" y**b, z**b", file=stderr)
print(y**b, file=stderr)
print(~(z**b), file=stderr)

print(" b**y, b**z", file=stderr)  # wolfram d^k/dx^k 4.0^x where x = 3.0
print(b**y, file=stderr)
print(~(b**z), file=stderr)

print(" (y + 1)**(y - 1), (z + 1)**(z - 1)", file=stderr)  # wolfram: d^k/dx^k (x + 1)^(x - 1) where x = 3.0
print((y + 1)**(y - 1), file=stderr)
print(~((z + 1)**(z - 1)), file=stderr)

print(" sqrt", file=stderr)
print(Dual.get(b).var.sqrt, file=stderr)
print(~Series.get(order, b).var.sqrt, file=stderr)

print(" exp", file=stderr)
print(Dual.get(b).var.exp, file=stderr)
print(~Series.get(order, b).var.exp, file=stderr)

print(" ln", file=stderr)
print(Dual.get(b).var.ln, file=stderr)
print(~Series.get(order, b).var.ln, file=stderr)

si = Dual.get(pi / a).var.sin
co = Dual.get(pi / a).var.cos
sine, cosine = Series.get(order, pi / a).var.sin_cos

print(f" sin(pi / {a:.1})", file=stderr)
print(si, file=stderr)
print(~sine, file=stderr)

print(f" cos(pi / {a:.1})", file=stderr)
print(co, file=stderr)
print(~cosine, file=stderr)

print(f" tan(pi / {b:.1})", file=stderr)
print(Dual.get(pi / b).var.tan, file=stderr)
print(~Series.get(order, pi / b).var.tan, file=stderr)

si = Dual.get(pi / a).var.sinh
co = Dual.get(pi / a).var.cosh
sine, cosine = Series.get(order, pi / a).var.sinh_cosh

print(f" sinh(pi / {a:.1})", file=stderr)
print(si, file=stderr)
print(~sine, file=stderr)

print(f" cosh(pi / {a:.1})", file=stderr)
print(co, file=stderr)
print(~cosine, file=stderr)

print(f" tanh(pi / {b:.1})", file=stderr)
print(Dual.get(pi / b).var.tanh, file=stderr)
print(~Series.get(order, pi / b).var.tanh, file=stderr)

# print(f" sec(pi / {b:.1})^2", file=stderr)
# print(~Series.get(order, pi / b).var.sec2, file=stderr)

print(f" asin({c:.1})", file=stderr)
print(Dual.get(c).var.asin, file=stderr)
print(~Series.get(order, c).var.asin, file=stderr)

print(f" acos({c:.1})", file=stderr)
print(Dual.get(c).var.acos, file=stderr)
print(~Series.get(order, c).var.acos, file=stderr)

print(f" atan({c:.1})", file=stderr)
print(Dual.get(c).var.atan, file=stderr)
print(~Series.get(order, c).var.atan, file=stderr)

print(f" asinh({c:.1})", file=stderr)
print(Dual.get(c).var.asinh, file=stderr)
print(~Series.get(order, c).var.asinh, file=stderr)

print(f" acosh({a:.1})", file=stderr)
print(Dual.get(a).var.acosh, file=stderr)
print(~Series.get(order, a).var.acosh, file=stderr)

print(f" atanh({c:.1})", file=stderr)
print(Dual.get(c).var.atanh, file=stderr)
print(~Series.get(order, c).var.atanh, file=stderr)

print("", file=stderr)

def analyze(model, x0, x1, steps):
    for k in range(steps):
        w_x = Dual.get(x_step(x0, x1, steps, k)).var
        w_f = model(w_x)
        print(f"{w_x.val:.6e} {w_f}")
        yield w_x.val, w_f.val, w_f.der

lower = float(argv[2])
upper = float(argv[3])
assert upper > lower
n_steps = int(argv[4])
assert n_steps > 0

for result in analyze(composite1, lower, upper, steps=n_steps):
    pass
