#!/usr/bin/env python3

from sys import argv, stderr
from math import acos
from ad import Series, Dual
from functions import playground

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

print("b**2", file=stderr)
print(Dual.get(b).var**2, file=stderr)
print(~Series.get(order, b).var**2, file=stderr)

print("b**-2", file=stderr)
print(Dual.get(b).var**-2, file=stderr)
print(~Series.get(order, b).var**-2, file=stderr)

print("b**0.5", file=stderr)
print(Dual.get(b).var**c, file=stderr)
print(~Series.get(order, b).var**c, file=stderr)

print("b**-0.5", file=stderr)
print(Dual.get(b).var**-c, file=stderr)
print(~Series.get(order, b).var**-c, file=stderr)

print(f" exp({b:.1f})", file=stderr)
print(Dual.get(b).var.exp, file=stderr)
print(~Series.get(order, b).var.exp, file=stderr)

print(f" ln({b:.1f})", file=stderr)
print(Dual.get(b).var.ln, file=stderr)
print(~Series.get(order, b).var.ln, file=stderr)

print(f" sin(π / {a:.0f})", file=stderr)
print(Dual.get(pi / a).var.sin, file=stderr)
print(~Series.get(order, pi / a).var.sin, file=stderr)

print(f" cos(π / {a:.0f})", file=stderr)
print(Dual.get(pi / a).var.cos, file=stderr)
print(~Series.get(order, pi / a).var.cos, file=stderr)

print(f" tan(π / {b:.0f})", file=stderr)
print(Dual.get(pi / b).var.tan, file=stderr)
print(~Series.get(order, pi / b).var.tan, file=stderr)

print(f" sinh(π / {a:.0f})", file=stderr)
print(Dual.get(pi / a).var.sinh, file=stderr)
print(~Series.get(order, pi / a).var.sinh, file=stderr)

print(f" cosh(π / {a:.0f})", file=stderr)
print(Dual.get(pi / a).var.cosh, file=stderr)
print(~Series.get(order, pi / a).var.cosh, file=stderr)

print(f" tanh(π / {b:.0f})", file=stderr)
print(Dual.get(pi / b).var.tanh, file=stderr)
print(~Series.get(order, pi / b).var.tanh, file=stderr)

print("", file=stderr)

def analyze(model, x0, x1, steps):
    for k in range(steps):
        w_x = Dual.get(x0 + k * (x1 - x0) / (steps - 1)).var
        w_f = model(w_x)
        print(f"{w_x.val:.6e} {w_f}")
        yield w_x.val, w_f.val, w_f.der

lower = float(argv[2])
upper = float(argv[3])
assert upper > lower
n_steps = int(argv[4])
assert n_steps > 0

for result in analyze(playground, lower, upper, steps=n_steps):
    pass
