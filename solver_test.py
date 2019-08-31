
#  Plotting
#  ./models.py BI 9 -5 5 1001 13 1e-12 1e-12 | ./plotMany.py 5 10 >/dev/null
#  Unit Testing
#  pytest --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py -v
from playground import bisect, newton, analyze, Solver, Mode
import pytest

order = 6
ε = 1.0e-12  # small error
δ = 1.0e-6  # small amount
xa = 3.0
xb = 4.0
xc = (xa + xb) / 2.0
f_tol = ε
x_tol = ε
max_it = 101
model = lambda x: x**3 - 4 * x**2 + 3 * x - 2
n_points = 101
plot_max = 8

a_max, b_max, x_max = 0.0, 1.0, +4.514162296e-01
a_infl, b_infl, x_infl = 1.0, 2.0, +1.333333333e+00
a_min, b_min, x_min = 2.0, 3.0, +2.215250437e+00
a_root, b_root, x_root = 3.0, 4.0, +3.269530842e+00

def test_bail():
    max_iterations = 1
    assert bisect(model, a_root, b_root, εf=f_tol, εx=x_tol, limit=max_iterations).count == max_iterations
    assert newton(model, (a_root + b_root) / 2.0, εf=f_tol, εx=x_tol, limit=max_iterations).count == max_iterations

@pytest.mark.parametrize("a, b, mode, target_f",
                         [(a_max, b_max, Mode.MIN_MAX, x_max),
                          (a_infl, b_infl, Mode.INFLECT, x_infl),
                          (a_min, b_min, Mode.MIN_MAX, x_min),
                          (a_root, b_root, Mode.ROOT___, x_root)])
def test_bisect(a, b, mode, target_f):
    result = bisect(model, a, b, εf=f_tol, εx=x_tol, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - result.x) < δ

@pytest.mark.parametrize("a, b, mode, target_f",
                         [(a_max, b_max, Mode.MIN_MAX, x_max),
                          (a_infl, b_infl, Mode.INFLECT, x_infl),
                          (a_min, b_min, Mode.MIN_MAX, x_min),
                          (a_root, b_root, Mode.ROOT___, x_root)])
def test_newton(a, b, mode, target_f):
    result = newton(model, (a + b) / 2.0, εf=f_tol, εx=x_tol, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < δ
    assert abs(result.f) < ε
    assert abs(target_f - result.x) < δ

def test_analysis_na():
    results = []
    for result in analyze(model, Solver.NA, - plot_max, plot_max, n_points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 0

@pytest.mark.parametrize("solver", [Solver.BI, Solver.NT])
def test_analysis(solver):
    results = []
    for result in analyze(model, solver, - plot_max, plot_max, n_points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 4
    assert results[0].mode == Mode.MIN_MAX.name
    assert abs(results[0].x - x_max) < δ
    assert results[1].mode == Mode.INFLECT.name
    assert abs(results[1].x - x_infl) < δ
    assert results[2].mode == Mode.MIN_MAX.name
    assert abs(results[2].x - x_min) < δ
    assert results[3].mode == Mode.ROOT___.name
    assert abs(results[3].x - x_root) < δ
