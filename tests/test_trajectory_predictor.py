import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.trajectory_predictor import TrajectoryPredictor


def helix_point(t, r=5.0, c=1.0):
    return (r * math.cos(t), r * math.sin(t), c * t)


def helix_deriv(t, r=5.0, c=1.0):
    return (-r * math.sin(t), r * math.cos(t), c)


def helix_deriv2(t, r=5.0, c=1.0):
    return (-r * math.cos(t), -r * math.sin(t), 0.0)


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def test_predicts_helix_position_accurately():
    r, c = 5.0, 1.0
    kappa_true = r / (r * r + c * c)
    tau_true = c / (r * r + c * c)

    t0 = 2.0
    P0 = helix_point(t0, r, c)
    v0 = helix_deriv(t0, r, c)
    a0 = helix_deriv2(t0, r, c)

    predictor = TrajectoryPredictor(substeps_per_step=50)
    dt = 0.5
    preds = predictor.predict(P0, v0, a0, kappa_true, tau_true, dt=dt, steps=4)

    for i, pred in enumerate(preds, start=1):
        t_true = t0 + i * dt
        true_pos = helix_point(t_true, r, c)
        err = dist(pred, true_pos)
        assert err < 1e-4, f"krok {i}: blad {err} za duzy"


def test_curvature_aware_prediction_beats_linear_extrapolation_on_curve():
    r, c = 5.0, 1.0
    kappa_true = r / (r * r + c * c)
    tau_true = c / (r * r + c * c)
    t0 = 1.0
    P0 = helix_point(t0, r, c)
    v0 = helix_deriv(t0, r, c)
    a0 = helix_deriv2(t0, r, c)

    predictor = TrajectoryPredictor(substeps_per_step=50)
    dt = 0.3
    steps = 5
    preds = predictor.predict(P0, v0, a0, kappa_true, tau_true, dt=dt, steps=steps)

    speed = math.sqrt(sum(x * x for x in v0))
    direction = tuple(x / speed for x in v0)
    linear_preds = [
        tuple(P0[j] + direction[j] * speed * dt * k for j in range(3))
        for k in range(1, steps + 1)
    ]

    curv_errs, lin_errs = [], []
    for i in range(steps):
        t_true = t0 + (i + 1) * dt
        true_pos = helix_point(t_true, r, c)
        curv_errs.append(dist(preds[i], true_pos))
        lin_errs.append(dist(linear_preds[i], true_pos))

    # na zakrzywionym torze predykcja oparta o krzywizne MUSI bic zwykla
    # ekstrapolacje liniowa, i to wyraznie, zwlaszcza na dalszym horyzoncie
    assert curv_errs[-1] < lin_errs[-1] / 5


def test_straight_line_falls_back_to_linear_extrapolation():
    predictor = TrajectoryPredictor()
    P0 = (0.0, 0.0, 0.0)
    v0 = (2.0, 0.0, 0.0)
    a0 = (0.0, 0.0, 0.0)  # brak przyspieszenia -- linia prosta
    preds = predictor.predict(P0, v0, a0, curvature=0.0, torsion=0.0, dt=1.0, steps=3)
    expected = [(2.0, 0.0, 0.0), (4.0, 0.0, 0.0), (6.0, 0.0, 0.0)]
    for p, e in zip(preds, expected):
        assert dist(p, e) < 1e-9


def test_zero_velocity_returns_current_position():
    predictor = TrajectoryPredictor()
    P0 = (1.0, 2.0, 3.0)
    preds = predictor.predict(P0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.1, 0.1, dt=1.0, steps=3)
    for p in preds:
        assert p == P0


def test_rejects_invalid_substeps():
    try:
        TrajectoryPredictor(substeps_per_step=0)
        assert False, "powinno rzucic ValueError"
    except ValueError:
        pass
