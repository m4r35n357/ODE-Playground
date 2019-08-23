
#  pytest -v --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py
from ad import Dual
from playground import bisect, falsi, secant, newton, householder, analyze, Solver

# Context.places = 6
order = 6
ε = 1.0e-12  # small error
δ = 1.0e-6  # small amount
xa = 1.0
xb = -1.0
xc = (xa + xb) / 2.0
f_tol = ε
x_tol = ε
max_it = 101
target_f = 0.0
model = lambda x: (x - 1)**2 / (x.cosh + 1).ln - 1
points = 101

def test_bisect():
    result = bisect(model, xa, xb, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_false_position():
    result = falsi(model, xa, xb, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, illinois=False, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_false_position_illinois():
    result = falsi(model, xa, xb, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_secant():
    result = secant(model, xa, xb, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_newton():
    result = newton(model, xc, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_householder_1():
    result = householder(model, xc, 2, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_householder_2():
    result = householder(model, xc, 3, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_householder_3():
    result = householder(model, xc, 4, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_householder_4():
    result = householder(model, xc, 5, y=target_f, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - model(Dual.get(result.x)).val) < ε

def test_analysis_bisection():
    results = []
    for result in analyze(model, Solver.BI, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_false_position():
    results = []
    for result in analyze(model, Solver.FP, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_false_position_illinois():
    results = []
    for result in analyze(model, Solver.FI, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_secant():
    results = []
    for result in analyze(model, Solver.SC, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_newton():
    results = []
    for result in analyze(model, Solver.NT, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_householder_newton():
    results = []
    for result in analyze(model, Solver.H1, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_householder_halley():
    results = []
    for result in analyze(model, Solver.H2, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_householder_3():
    results = []
    for result in analyze(model, Solver.H3, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6

def test_analysis_householder_4():
    results = []
    for result in analyze(model, Solver.H4, -8.0, 8.0, points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6
