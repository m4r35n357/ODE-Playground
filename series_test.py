#!/usr/bin/env python3

from sys import argv, stderr
from math import pi, e
from ad import Series, Dual
from functions import playground, composite1, composite2, septic, lorentz

order = int(argv[1])
assert order > 1

a = 3.0
b = 4.0
c = 0.5

w, x = Dual.get(b).var, Series.get(order, b).var
y, z = Dual.get(a).var, Series.get(order, a).var

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

print(" 1 / y, 1 / z", file=stderr)
print((1 / y), file=stderr)
print(~(1 / z), file=stderr)

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
print(w**2, file=stderr)
print(~(x**2), file=stderr)

print("b**-2", file=stderr)
print(w**-2, file=stderr)
print(~(x**-2), file=stderr)

print("b**0.5", file=stderr)
print(w**c, file=stderr)
print(~(x**c), file=stderr)

print("b**-0.5", file=stderr)
print(w**-c, file=stderr)
print(~(x**-c), file=stderr)

print("b**0", file=stderr)
print(w**0, file=stderr)
print(~(x**0), file=stderr)

print("b**0.0", file=stderr)
print(w**0.0, file=stderr)
print(~(x**0.0), file=stderr)

print(f" exp({b:.1f})", file=stderr)
print(w.exp, file=stderr)
print(~(x.exp), file=stderr)

print(f" ln({b:.1f})", file=stderr)
print(w.ln, file=stderr)
print(~(x.ln), file=stderr)

w, x = Dual.get(pi / b).var, Series.get(order, pi / b).var
y, z = Dual.get(pi / a).var, Series.get(order, pi / a).var

print(f" sin(π/{a:.0f})", file=stderr)
print(y.sin, file=stderr)
print(~(z.sin), file=stderr)

print(f" cos(π/{a:.0f})", file=stderr)
print(y.cos, file=stderr)
print(~(z.cos), file=stderr)

print(f" tan(π/{b:.0f})", file=stderr)
print(w.tan, file=stderr)
print(~(x.tan), file=stderr)

print(f" sinh(π/{a:.0f})", file=stderr)
print(y.sinh, file=stderr)
print(~(z.sinh), file=stderr)

print(f" cosh(π/{a:.0f})", file=stderr)
print(y.cosh, file=stderr)
print(~(z.cosh), file=stderr)

print(f" tanh(π/{b:.0f})", file=stderr)
print(w.tanh, file=stderr)
print(~(x.tanh), file=stderr)

print("", file=stderr)

#  Zero identities
w, x = Dual.get(b).var, Series.get(order, b).var
y, z = Dual.get(a).var, Series.get(order, a).var

print("1.0 / b**2 - b**-2.0", file=stderr)
print(1.0 / w**2 - w**-2.0, file=stderr)
print(~(1.0 / x**2 - x**-2.0), file=stderr)

print("(b**2)**0.5 - |b|", file=stderr)
print((w**2)**0.5 - abs(w), file=stderr)
print(~((x**2)**0.5) - abs(x), file=stderr)

print("exp(b) - e**b", file=stderr)
print((w.exp - e**w), file=stderr)
print((x.exp - e**x), file=stderr)

print("ln(e**b) - b", file=stderr)
print((w.exp.ln - w), file=stderr)
print((x.exp.ln - x), file=stderr)

print(f" sin(3.0 * a) - 3.0 * sin(a) + 4.0 * sin(a)**3", file=stderr)
print((3 * y).sin - 3.0 * y.sin + 4.0 * y.sin**3, file=stderr)
print(~((3 * z).sin - 3.0 * z.sin + 4.0 * z.sin**3), file=stderr)

print(f" cos(3.0 * a) + 3.0 * cos(a) - 4.0 * cos(a)**3", file=stderr)
print((3 * y).cos + 3.0 * y.cos - 4.0 * y.cos**3, file=stderr)
print(~((3 * z).cos + 3.0 * z.cos - 4.0 * z.cos**3), file=stderr)

print(f" tan(b) - sin(b) / cos(b)", file=stderr)
print(w.tan - w.sin / w.cos, file=stderr)
print(~(x.tan - x.sin / x.cos), file=stderr)

print(f" tanh(b) - sinh(b) / cosh(b)", file=stderr)
print(w.tanh - w.sinh / w.cosh, file=stderr)
print(~(x.tanh - x.sinh / x.cosh), file=stderr)

print(f" exp vs sinh", file=stderr)
print(0.5 * (y.exp - 1.0 / y.exp) - y.sinh, file=stderr)
print(~(0.5 * (z.exp - 1.0 / z.exp) - z.sinh), file=stderr)

print(f" exp vs cosh", file=stderr)
print(0.5 * (y.exp + 1.0 / y.exp) - y.cosh, file=stderr)
print(~(0.5 * (z.exp + 1.0 / z.exp) - z.cosh), file=stderr)

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

# For Python console:
'''
from ad import *

a = Series.get(5, 3.0)
b = Series.get(5, 5.0)
c = Series.get(5, 7.0)
x = Series.get(5, 2.0)

y = a * x**3 - b * x**2 + c * x - 5
print(~y)

y = a * x.var**3 - b * x.var**2 + c * x.var - 5
print(~y)

y = a.var * x**3 - b * x**2 + c * x - 5
print(~y)

y = a * x**3 - b.var * x**2 + c * x - 5
print(~y)

y = a * x**3 - b * x**2 + c.var * x - 5
print(~y)
'''
