#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Plotting
#  ./models.py BI 9 -5 5 1001 13 1e-12 1e-12 | ./plotMany.py 5 10 >/dev/null
#  Unit Testing
#  pytest --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py -v
from math import sqrt
from playground import bisect, newton, analyze, Solver, Mode, Sense
import pytest

order = 6
ε = 1.0e-9  # small error
δ = 1.0e-6  # small amount
xa, xb = 3.0, 4.0
xc = (xa + xb) / 2.0
f_tol = x_tol = ε
max_it = 101
n_points = 1000
plot_min, plot_max = - 8.0, 8.0

model = lambda x: x ** 3 - 4 * x ** 2 + 3 * x - 2
a_max, b_max, x_max = 0.0, 1.0, +4.514162296e-01
a_infl, b_infl, x_infl = 1.0, 2.0, +1.333333333e+00
a_min, b_min, x_min = 2.0, 3.0, +2.215250437e+00
a_root, b_root, x_root = 3.0, 4.0, +3.269530842e+00

def test_bail():
    max_iterations = 1
    assert bisect(model, a_root, b_root, εf=f_tol, εx=x_tol, limit=max_iterations).count == max_iterations
    assert newton(model, (a_root + b_root) / 2.0, εf=f_tol, εx=x_tol, limit=max_iterations).count == max_iterations

def test_basic_solve():
    square = lambda x: x**2 - a
    a = 2.0
    result = bisect(square, 1.0, 2.0, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert abs(result.x - sqrt(a)) < ε
    result = newton(square, 1.0, εf=f_tol, εx=x_tol, limit=max_it, debug=True)
    assert abs(result.x - sqrt(a)) < ε

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
    for result in analyze(model, Solver.NA, plot_min, plot_max, n_points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 0

@pytest.mark.parametrize("solver", [Solver.BI, Solver.NT])
def test_analysis_simple(solver):
    results = []
    for result in analyze(model, solver, plot_min, plot_max, n_points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 4
    assert (results[0].mode == Mode.MIN_MAX.name) and (results[0].sense == Sense.DECREASING.value)
    assert abs(results[0].x - x_max) < δ
    assert (results[1].mode == Mode.INFLECT.name) and (results[1].sense == Sense.INCREASING.value)
    assert abs(results[1].x - x_infl) < δ
    assert (results[2].mode == Mode.MIN_MAX.name) and (results[2].sense == Sense.INCREASING.value)
    assert abs(results[2].x - x_min) < δ
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.INCREASING.value)
    assert abs(results[3].x - x_root) < δ

@pytest.mark.parametrize("solver", [Solver.BI, Solver.NT])
def test_analysis_septic(solver):
    results = []
    for result in analyze(lambda a: (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6),
                          solver, plot_min, plot_max, n_points, f_tol, x_tol, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 18
    # NT  x: -7.000e+00  δx: -0.000e+00  f: +0.000e+00  / ROOT___ 4
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.INCREASING.value)
    # NT  x: -6.321e+00  δx: +3.286e-16  f: +1.933e-11  \ MIN_MAX 4
    assert (results[1].mode == Mode.MIN_MAX.name) and (results[1].sense == Sense.DECREASING.value)
    # NT  x: -5.514e+00  δx: +9.065e-17  f: -3.638e-12  / INFLECT 4
    assert (results[2].mode == Mode.INFLECT.name) and (results[2].sense == Sense.INCREASING.value)
    # NT  x: -5.000e+00  δx: +0.000e+00  f: +0.000e+00  \ ROOT___ 4
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.DECREASING.value)
    # NT  x: -3.901e+00  δx: +6.086e-14  f: -6.934e-10  / MIN_MAX 3
    assert (results[4].mode == Mode.MIN_MAX.name) and (results[4].sense == Sense.INCREASING.value)
    # NT  x: -2.836e+00  δx: +5.234e-17  f: +4.547e-13  \ INFLECT 4
    assert (results[5].mode == Mode.INFLECT.name) and (results[5].sense == Sense.DECREASING.value)
    # NT  x: -2.000e+00  δx: -0.000e+00  f: +0.000e+00  / ROOT___ 4
    assert (results[6].mode == Mode.ROOT___.name) and (results[6].sense == Sense.INCREASING.value)
    # NT  x: -1.192e+00  δx: +1.542e-16  f: +5.684e-13  \ MIN_MAX 4
    assert (results[7].mode == Mode.MIN_MAX.name) and (results[7].sense == Sense.DECREASING.value)
    # NT  x: -3.756e-01  δx: -1.182e-17  f: +5.684e-14  / INFLECT 4
    assert (results[8].mode == Mode.INFLECT.name) and (results[8].sense == Sense.INCREASING.value)
    # NT  x: +0.000e+00  δx: +9.808e-19  f: +1.236e-15  \ ROOT___ 4
    assert (results[9].mode == Mode.ROOT___.name) and (results[9].sense == Sense.DECREASING.value)
    # NT  x: +5.159e-01  δx: -3.326e-17  f: +9.948e-14  / MIN_MAX 4
    assert (results[10].mode == Mode.MIN_MAX.name) and (results[10].sense == Sense.INCREASING.value)
    # NT  x: +1.000e+00  δx: -0.000e+00  f: +0.000e+00  / ROOT___ 4
    assert (results[11].mode == Mode.ROOT___.name) and (results[11].sense == Sense.INCREASING.value)
    # NT  x: +1.557e+00  δx: +3.274e-17  f: +2.274e-13  \ INFLECT 4
    assert (results[12].mode == Mode.INFLECT.name) and (results[12].sense == Sense.DECREASING.value)
    # NT  x: +2.294e+00  δx: +1.363e-16  f: +9.095e-13  \ MIN_MAX 4
    assert (results[13].mode == Mode.MIN_MAX.name) and (results[13].sense == Sense.DECREASING.value)
    # NT  x: +3.000e+00  δx: +0.000e+00  f: +0.000e+00  \ ROOT___ 4
    assert (results[14].mode == Mode.ROOT___.name) and (results[14].sense == Sense.DECREASING.value)
    # NT  x: +4.311e+00  δx: -1.912e-16  f: +7.276e-12  / INFLECT 4
    assert (results[15].mode == Mode.INFLECT.name) and (results[15].sense == Sense.INCREASING.value)
    # NT  x: +5.175e+00  δx: +1.166e-16  f: -7.276e-12  / MIN_MAX 4
    assert (results[16].mode == Mode.MIN_MAX.name) and (results[16].sense == Sense.INCREASING.value)
    # NT  x: +6.000e+00  δx: -8.882e-16  f: +9.145e-11  / ROOT___ 4
    assert (results[17].mode == Mode.ROOT___.name) and (results[17].sense == Sense.INCREASING.value)

