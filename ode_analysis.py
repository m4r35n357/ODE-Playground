import datetime
from sys import stderr
from collections import namedtuple
from matplotlib import pyplot
from dual import Context, Dual, Result, Solver, Sense

class Matrix3x3(namedtuple('Matrix3x3Type', ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'])):
    def __str__(self):
        return f'{self.a:+.{Context.places}e}, {self.b:+.{Context.places}e}, {self.c:+.{Context.places}e}\n' \
               f'{self.d:+.{Context.places}e}, {self.e:+.{Context.places}e}, {self.f:+.{Context.places}e}\n' \
               f'{self.g:+.{Context.places}e}, {self.h:+.{Context.places}e}, {self.i:+.{Context.places}e}'

def function_3(model_a, model_b, model_c, x, y, z):
    a, b, c = Dual.get(x), Dual.get(y), Dual.get(z)
    return model_a(a, b, c).val, model_b(a, b, c).val, model_c(a, b, c).val

def jacobian_3x3(model_a, model_b, model_c, x, y, z):
    a, b, c = Dual.get(x), Dual.get(y), Dual.get(z)
    return Matrix3x3(a=model_a(a.var, b, c).der, b=model_a(a, b.var, c).der, c=model_a(a, b, c.var).der,
                     d=model_b(a.var, b, c).der, e=model_b(a, b.var, c).der, f=model_b(a, b, c.var).der,
                     g=model_c(a.var, b, c).der, h=model_c(a, b.var, c).der, i=model_c(a, b, c.var).der)

def determinant_3x3(m):
    cofactors = Matrix3x3(a=(m.e * m.i - m.f * m.h), b=-(m.d * m.i - m.f * m.g), c=(m.d * m.h - m.e * m.g),
                          d=-(m.b * m.i - m.c * m.h), e=(m.a * m.i - m.c * m.g), f=-(m.a * m.h - m.b * m.g),
                          g=(m.b * m.f - m.c * m.e), h=-(m.a * m.f - m.c * m.d), i=(m.a * m.e - m.b * m.d))
    determinant =  m.a * cofactors.a + m.b * cofactors.b + m.c * cofactors.c
    return determinant, cofactors

def invert_3x3(m):
    d, c = determinant_3x3(m)
    assert d != 0.0, f"determinant = {d}"
    return Matrix3x3(a=c.a/d, b=c.d/d, c=c.g/d,
                     d=c.b/d, e=c.e/d, f=c.h/d,
                     g=c.c/d, h=c.f/d, i=c.i/d)

def equilibrium(model_a, model_b, model_c, x, y, z, εf=1e-12, max_it=100):
    fx, fy, fz = function_3(model_a, model_b, model_c, x, y, z)
    count = 0
    print(f'{count}, {x:+.{Context.places}e}, {y:+.{Context.places}e}, {z:+.{Context.places}e}, '
          f'{fx:+.{Context.places}e}, {fy:+.{Context.places}e}, {fz:+.{Context.places}e}')
    while abs(fx) + abs(fy) + abs(fz) > εf and count < max_it:
        j_1 = invert_3x3(jacobian_3x3(model_a, model_b, model_c, x, y, z))
        x -= j_1.a * fx + j_1.b * fy + j_1.c * fz
        y -= j_1.d * fx + j_1.e * fy + j_1.f * fz
        z -= j_1.g * fx + j_1.h * fy + j_1.i * fz
        fx, fy, fz = function_3(model_a, model_b, model_c, x, y, z)
        count += 1
        print(f'{count}, {x:+.{Context.places}e}, {y:+.{Context.places}e}, {z:+.{Context.places}e}, '
              f'{fx:+.{Context.places}e}, {fy:+.{Context.places}e}, {fz:+.{Context.places}e}')
    return x, y, z, fx, fy, fz

def characteristic_equation(m, λ):
    return determinant_3x3(Matrix3x3(a=m.a - λ, b=m.b, c=m.c,
                                     d=m.d, e=m.e - λ, f=m.f,
                                     g=m.g, h=m.h, i=m.i - λ))[0]

def bisect_ce(j, λa, λb, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    a, b, c = Dual.get(λa), Dual.get(λb), Dual.get()
    ce_sign = characteristic_equation(j, Dual.get(λa)).val
    δx = ce = i = 1
    while i <= limit and abs(ce) > εf or abs(δx) > εx:
        c = 0.5 * (a + b)
        ce = characteristic_equation(j, c).val
        if ce_sign * ce < 0.0:
            b = c
        elif ce_sign * ce > 0.0:
            a = c
        else:
            break
        δx = b.val - a.val
        i += 1
        if debug:
            print(Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=ce, δx=δx), file=stderr)
    return Result(method=Solver.BI.name, count=i-1, sense=sense.value, x=c.val, f=ce, δx=δx)

def newton_ce(j, λ0, εf=1e-12, εx=1e-12, limit=101, sense=Sense.FLAT, debug=False):
    λ, ce = Dual.get(λ0).var, Dual.get(1)
    δλ = i = 1
    while i <= limit and abs(ce.val) > εf or abs(δλ) > εx:
        ce = characteristic_equation(j, λ)
        δλ = - ce.val / ce.der
        λ += δλ
        i += 1
        if debug:
            print(Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=λ.val, f=ce.val, δx=δλ), file=stderr)
    return Result(method=Solver.NT.name, count=i-1, sense=sense.value, x=λ.val, f=ce.val, δx=δλ)

def analyze_ce(j, method, λa, λb, steps, εf=1e-12, εx=1e-12, limit=101, debug=False):
    λ_prev = ce_prev = None
    step = (λb - λa) / (steps - 1)
    for k in range(steps):
        λ = Dual.get(λa + k * step).var
        ce = characteristic_equation(j, λ)
        if k > 0:
            if ce_prev * ce.val < 0.0:
                sense = Sense.DECREASING if ce_prev > ce.val else Sense.INCREASING
                if method == Solver.BI:
                    yield bisect_ce(j, λ.val, λ_prev, εf, εx, limit=limit, sense=sense, debug=debug)
                if method == Solver.NT:
                    yield newton_ce(j, λ.val, εf, εx, limit, sense=sense, debug=debug)
        λ_prev = λ.val
        ce_prev = ce.val

def _plot(model_a, model_b, model_c, x, y, z, λ_min, λ_max, steps, ce_min, ce_max, εf, εx, limit, nt, debug):
    solver = Solver.NT if nt else Solver.BI
    j = jacobian_3x3(model_a, model_b, model_c, x, y, z)
    print(f'Jacobian({x}, {y}, {z})\n{j}')
    for result in analyze_ce(j, solver, λ_min, λ_max, steps, εf, εx, limit, debug):
        # noinspection PyTypeChecker
        if result.count < 101:
            print(f'eigenvalue: {result.x}\n{result}')
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel('λ', color='.2')
    ax1.set_ylabel(f'Characteristic Equation and first derivative', color='.2')
    ax1.set_xlim(λ_min, λ_max)
    ax1.set_ylim(ce_min, ce_max)
    colour = ['k-', 'r-', 'g-', 'b-', 'y-', 'c-', 'm-', 'r:', 'g:', 'b:', 'y:', 'c:', 'm:']
    data = [[] for _ in range(3)]
    step = (λ_max - λ_min) / (steps - 1)
    for k in range(steps):
        λ = Dual.get(λ_min + k * step).var
        ce = characteristic_equation(j, λ)
        for d_term, p_term in zip(data, [λ.val] + [ce.val] + [ce.der]):
            d_term.append(p_term)
    for c in reversed(range(1, 3)):
        ax1.plot(data[0], data[c], f'{colour[c - 1]}', linewidth=2 if c == 1 else 1, markersize=0, label=c-1)
    ax1.legend(loc='lower right')

def plot_lambda(model_a, model_b, model_c, x, y, z, λ_min=-50.0, λ_max=50.0, steps=1000, ce_min=-5000.0, ce_max=5000.0,
                εf=1e-12, εx=1e-12, limit=101, nt=True, debug=False):
    _plot(model_a, model_b, model_c, x, y, z, λ_min=λ_min, λ_max=λ_max, steps=steps, ce_min=ce_min, ce_max=ce_max,
                εf=εf, εx=εx, limit=limit, nt=nt, debug=debug)
    pyplot.show()

def save_lambda(model_a, model_b, model_c, x, y, z, λ_min=-50.0, λ_max=50.0, steps=1000, ce_min=-1000.0, ce_max=1000.0,
                εf=1e-12, εx=1e-12, limit=101, nt=True, debug=False):
    _plot(model_a, model_b, model_c, x, y, z, λ_min=λ_min, λ_max=λ_max, steps=steps, ce_min=-ce_min, ce_max=ce_max,
                εf=εf, εx=εx, limit=limit, nt=nt, debug=debug)
    pyplot.savefig(f'/tmp/plot-{datetime.datetime.now()}.png')
