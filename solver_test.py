#
#  (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov=playground --cov-report html:cov_html ad_test.py solver_test.py -v
from math import sqrt
from playground import bisect, newton, analyze, Solver, Mode, Sense
import pytest

order = 6
εf = 1.0e-9
εx = 1.0e-9
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
    assert bisect(lambda a: a**2 - 2, 1.0, 2.0, εf=εf, εx=εx, limit=max_iterations).count == max_iterations
    assert newton(lambda a: a**2 - 2, 1.5, εf=εf, εx=εx, limit=max_iterations).count == max_iterations

def test_basic_solve():
    target = 2.0
    result = bisect(lambda a: a**2 - target, 1.0, 2.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx
    result = newton(lambda a: a**2 - target, 1.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx

@pytest.mark.parametrize('a, b, mode, target_x',
                         [(a_max, b_max, Mode.MIN_MAX, x_max),
                          (a_infl, b_infl, Mode.INFLECT, x_infl),
                          (a_min, b_min, Mode.MIN_MAX, x_min),
                          (a_root, b_root, Mode.ROOT___, x_root)])
def test_bisect(a, b, mode, target_x):
    result = bisect(model, a, b, εf=εf, εx=εx, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx

@pytest.mark.parametrize('a, b, mode, target_x',
                         [(a_max, b_max, Mode.MIN_MAX, x_max),
                          (a_infl, b_infl, Mode.INFLECT, x_infl),
                          (a_min, b_min, Mode.MIN_MAX, x_min),
                          (a_root, b_root, Mode.ROOT___, x_root)])
def test_newton(a, b, mode, target_x):
    result = newton(model, (a + b) / 2.0, εf=εf, εx=εx, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx

def test_analysis_na(capsys):
    captured = capsys.readouterr()
    results = []
    for result in analyze(model, Solver.NA, plot_min, plot_max, n_points, εf, εx, limit=max_it, order=order, console=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 0
    assert len(captured.out) == 0

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_cubic(solver):
    results = []
    for result in analyze(lambda a: a**3 + 3 * a**2 - 3,
                          solver, plot_min, plot_max, n_points, εf, εx, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6  # 3 roots, 1 maximum, 1 minimum, 1 inflection
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.MIN_MAX.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.ROOT___.name) and (results[2].sense == Sense.DECREASING.value)
    assert (results[3].mode == Mode.INFLECT.name) and (results[3].sense == Sense.INCREASING.value)
    assert (results[4].mode == Mode.MIN_MAX.name) and (results[4].sense == Sense.INCREASING.value)
    assert (results[5].mode == Mode.ROOT___.name) and (results[5].sense == Sense.INCREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_cos_cubic(solver):
    results = []
    for result in analyze(lambda a: a.cos - a**3,
                          solver, plot_min, plot_max, n_points, εf, εx, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 4
    assert (results[0].mode == Mode.MIN_MAX.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.INFLECT.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.MIN_MAX.name) and (results[2].sense == Sense.DECREASING.value)
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.DECREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_messy(solver):
    results = []
    for result in analyze(lambda a: (a - 1)**2 / (a.cosh + 1).ln - 1,
                          solver, plot_min, plot_max, n_points, εf, εx, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 6
    assert (results[0].mode == Mode.INFLECT.name) and (results[0].sense == Sense.DECREASING.value)
    assert (results[1].mode == Mode.INFLECT.name) and (results[1].sense == Sense.INCREASING.value)
    assert (results[2].mode == Mode.ROOT___.name) and (results[2].sense == Sense.DECREASING.value)
    assert (results[3].mode == Mode.MIN_MAX.name) and (results[3].sense == Sense.INCREASING.value)
    assert (results[4].mode == Mode.ROOT___.name) and (results[4].sense == Sense.INCREASING.value)
    assert (results[5].mode == Mode.INFLECT.name) and (results[5].sense == Sense.DECREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_septic(solver):
    results = []
    for result in analyze(lambda a: (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6),
                          solver, plot_min, plot_max, n_points, εf, εx, limit=max_it, order=order):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 18  # 7 roots, 3 maxima, 3 minima, 5 inflections
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.MIN_MAX.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.INFLECT.name) and (results[2].sense == Sense.INCREASING.value)
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.DECREASING.value)
    assert (results[4].mode == Mode.MIN_MAX.name) and (results[4].sense == Sense.INCREASING.value)
    assert (results[5].mode == Mode.INFLECT.name) and (results[5].sense == Sense.DECREASING.value)
    assert (results[6].mode == Mode.ROOT___.name) and (results[6].sense == Sense.INCREASING.value)
    assert (results[7].mode == Mode.MIN_MAX.name) and (results[7].sense == Sense.DECREASING.value)
    assert (results[8].mode == Mode.INFLECT.name) and (results[8].sense == Sense.INCREASING.value)
    assert (results[9].mode == Mode.ROOT___.name) and (results[9].sense == Sense.DECREASING.value)
    assert (results[10].mode == Mode.MIN_MAX.name) and (results[10].sense == Sense.INCREASING.value)
    assert (results[11].mode == Mode.ROOT___.name) and (results[11].sense == Sense.INCREASING.value)
    assert (results[12].mode == Mode.INFLECT.name) and (results[12].sense == Sense.DECREASING.value)
    assert (results[13].mode == Mode.MIN_MAX.name) and (results[13].sense == Sense.DECREASING.value)
    assert (results[14].mode == Mode.ROOT___.name) and (results[14].sense == Sense.DECREASING.value)
    assert (results[15].mode == Mode.INFLECT.name) and (results[15].sense == Sense.INCREASING.value)
    assert (results[16].mode == Mode.MIN_MAX.name) and (results[16].sense == Sense.INCREASING.value)
    assert (results[17].mode == Mode.ROOT___.name) and (results[17].sense == Sense.INCREASING.value)

