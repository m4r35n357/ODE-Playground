#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
#  Unit Testing
#  pytest --cov=ad --cov=ad --cov-report html:cov_html ad_test.py solver_test.py -v
from math import sqrt
from ad import Solver, Mode, Sense, s_bisect, s_newton, s_analyze, d_bisect, d_newton, d_analyze
import pytest

εf = 1.0e-9
εx = 1.0e-9
max_it = 101
n_points = 1000
plot_min, plot_max = - 8.0, 8.0

model = lambda x: x ** 3 - 4 * x ** 2 + 3 * x - 2

def test_bail():
    max_iterations = 1
    assert s_bisect(lambda a: a ** 2 - 2, 1.0, 2.0, εf=εf, εx=εx, limit=max_iterations).count == max_iterations
    assert s_newton(lambda a: a ** 2 - 2, 1.5, εf=εf, εx=εx, limit=max_iterations).count == max_iterations
    assert d_bisect(lambda a: a ** 2 - 2, 1.0, 2.0, εf=εf, εx=εx, limit=max_iterations).count == max_iterations
    assert d_newton(lambda a: a ** 2 - 2, 1.5, εf=εf, εx=εx, limit=max_iterations).count == max_iterations

def test_quadratic_solve_series():
    target = 2.0
    result = s_bisect(lambda a: a ** 2 - target, 1.0, 2.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx
    result = s_newton(lambda a: a ** 2 - target, 1.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx

def test_quadratic_solve_dual():
    target = 2.0
    result = d_bisect(lambda a: a ** 2 - target, 1.0, 2.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx
    result = d_newton(lambda a: a ** 2 - target, 1.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.x - sqrt(target)) < εx

@pytest.mark.parametrize('a, b, mode, target_x',
                         [(0.0, 1.0, Mode.MIN_MAX, +4.514162296e-01),
                          (1.0, 2.0, Mode.INFLECT, +1.333333333e+00),
                          (2.0, 3.0, Mode.MIN_MAX, +2.215250437e+00),
                          (3.0, 4.0, Mode.ROOT___, +3.269530842e+00)])
def test_cubic_solve_series(a, b, mode, target_x):
    result = s_bisect(model, a, b, εf=εf, εx=εx, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx
    result = s_newton(model, (a + b) / 2.0, εf=εf, εx=εx, limit=max_it, mode=mode, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx

@pytest.mark.parametrize('a, b, mode, target_x',
                         [(3.0, 4.0, Mode.ROOT___, +3.269530842e+00)])
def test_cubic_solve_dual(a, b, mode, target_x):
    result = d_bisect(model, a, b, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx
    result = d_newton(model, (a + b) / 2.0, εf=εf, εx=εx, limit=max_it, debug=True)
    assert result.count < max_it
    assert abs(result.δx) < εx
    assert abs(result.f) < εf
    assert abs(target_x - result.x) < εx

def test_analysis_na_series(capsys):
    captured = capsys.readouterr()
    results = []
    for result in s_analyze(model=model, method=Solver.NA, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, mode=Mode.ALL, console=False, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 0
    assert len(captured.out) == 0

def test_analysis_na_dual(capsys):
    captured = capsys.readouterr()
    results = []
    for result in d_analyze(model=model, method=Solver.NA, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, console=False, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 0
    assert len(captured.out) == 0

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_cubic_series(solver):
    results = []
    for result in s_analyze(model=lambda a: a ** 3 + 3 * a ** 2 - 3,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, mode=Mode.ALL, console=True, debug=False):
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
def test_analysis_cubic_dual(solver):
    results = []
    for result in d_analyze(model=lambda a: a ** 3 + 3 * a ** 2 - 3,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, console=True, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 3  # 3 roots
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.ROOT___.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.ROOT___.name) and (results[2].sense == Sense.INCREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_cos_cubic_series(solver):
    results = []
    for result in s_analyze(model=lambda a: a.cos - a ** 3,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, mode=Mode.ALL, console=True, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 4
    assert (results[0].mode == Mode.MIN_MAX.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.INFLECT.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.MIN_MAX.name) and (results[2].sense == Sense.DECREASING.value)
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.DECREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_cos_cubic_dual(solver):
    results = []
    for result in d_analyze(model=lambda a: a.cos - a ** 3,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, console=True, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 1
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.DECREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_messy_series(solver):
    results = []
    for result in s_analyze(model=lambda a: (a - 1) ** 2 / (a.cosh + 1).ln - 1,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, mode=Mode.ALL, console=True, debug=False):
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
def test_analysis_messy_dual(solver):
    results = []
    for result in d_analyze(model=lambda a: (a - 1) ** 2 / (a.cosh + 1).ln - 1,
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, console=True, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 2
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.DECREASING.value)
    assert (results[1].mode == Mode.ROOT___.name) and (results[1].sense == Sense.INCREASING.value)

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_septic_series(solver):
    results = []
    for result in s_analyze(model=lambda a: (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6),
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, mode=Mode.ALL, console=True, debug=False):
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

@pytest.mark.parametrize('solver', [Solver.BI, Solver.NT])
def test_analysis_septic_dual(solver):
    results = []
    for result in d_analyze(model=lambda a: (a + 7) * (5 + a) * (a + 2.0) * a * (1 - a) * (3.0 - a) * (a - 6),
                            method=solver, x0=plot_min, x1=plot_max, steps=n_points, εf=εf, εx=εx,
                            limit=max_it, console=True, debug=False):
        if result.count < max_it:
            results.append(result)
    assert len(results) == 7  # 7 roots
    assert (results[0].mode == Mode.ROOT___.name) and (results[0].sense == Sense.INCREASING.value)
    assert (results[1].mode == Mode.ROOT___.name) and (results[1].sense == Sense.DECREASING.value)
    assert (results[2].mode == Mode.ROOT___.name) and (results[2].sense == Sense.INCREASING.value)
    assert (results[3].mode == Mode.ROOT___.name) and (results[3].sense == Sense.DECREASING.value)
    assert (results[4].mode == Mode.ROOT___.name) and (results[4].sense == Sense.INCREASING.value)
    assert (results[5].mode == Mode.ROOT___.name) and (results[5].sense == Sense.DECREASING.value)
    assert (results[6].mode == Mode.ROOT___.name) and (results[6].sense == Sense.INCREASING.value)
