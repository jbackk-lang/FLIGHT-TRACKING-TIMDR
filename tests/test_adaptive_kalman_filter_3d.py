import csv
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.kalman_filter_3d import KalmanFilter3D
from core.adaptive_kalman_filter_3d import AdaptiveKalmanFilter3D

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _err(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _load_sharp_maneuver():
    path = os.path.join(DATA_DIR, "synthetic_maneuvers.csv")
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    rows = [r for r in reader if r["segment"] == "sharp_maneuver"]
    pts = [(float(r["x"]), float(r["y"]), float(r["z"])) for r in rows]
    true_pts = [(float(r["x_true"]), float(r["y_true"]), float(r["z_true"])) for r in rows]
    return pts, true_pts


def test_adaptive_q_reduces_error_during_maneuver_vs_static_kalman():
    pts, true_pts = _load_sharp_maneuver()
    static = KalmanFilter3D(process_var=0.01, measurement_var=5.0, initial_position=pts[0])
    adaptive = AdaptiveKalmanFilter3D(process_var=0.01, measurement_var=5.0, initial_position=pts[0])

    static_errs, adaptive_errs = [], []
    for p, tp in zip(pts, true_pts):
        se = static.update(p)
        ae = adaptive.update(p)
        static_errs.append(_err(se, tp))
        adaptive_errs.append(_err(ae, tp))

    # manewr jest w probkach ~40-50 (patrz generate_synthetic_maneuvers.py)
    maneuver_window = range(40, 55)
    static_mean = sum(static_errs[i] for i in maneuver_window) / len(maneuver_window)
    adaptive_mean = sum(adaptive_errs[i] for i in maneuver_window) / len(maneuver_window)
    assert adaptive_mean < static_mean * 0.7  # co najmniej 30% mniejszy blad w oknie manewru

    # na prostych odcinkach adaptacja NIE MOZE znaczaco pogorszyc sledzenia
    straight_window = list(range(0, 35)) + list(range(60, 80))
    static_straight = sum(static_errs[i] for i in straight_window) / len(straight_window)
    adaptive_straight = sum(adaptive_errs[i] for i in straight_window) / len(straight_window)
    assert adaptive_straight < static_straight * 1.2


def test_maneuvering_flag_activates_during_and_only_during_maneuver():
    pts, _ = _load_sharp_maneuver()
    adaptive = AdaptiveKalmanFilter3D(process_var=0.01, measurement_var=5.0, initial_position=pts[0])
    flags = []
    for p in pts:
        adaptive.update(p)
        flags.append(adaptive.last_maneuvering)

    # zaden alarm manewru na pierwszych 30 probkach czysto prostego lotu
    assert not any(flags[:30])
    # co najmniej jeden alarm gdzies w oknie manewru + krotki ogon po nim
    assert any(flags[40:55])


def test_measurement_quality_improves_tracking_during_noise_burst():
    random.seed(7)
    true_pts = [(2.0 * t, 0.0, 0.0) for t in range(60)]
    noisy_pts, quality = [], []
    for i, p in enumerate(true_pts):
        if 25 <= i < 35:
            noisy_pts.append(tuple(c + random.gauss(0, 5.0) for c in p))
            quality.append(0.15)
        else:
            noisy_pts.append(tuple(c + random.gauss(0, 0.3) for c in p))
            quality.append(1.0)

    ignore_quality = AdaptiveKalmanFilter3D(process_var=0.01, measurement_var=0.5, initial_position=noisy_pts[0])
    use_quality = AdaptiveKalmanFilter3D(process_var=0.01, measurement_var=0.5, initial_position=noisy_pts[0])

    e_ignore, e_use = [], []
    for p, tp, q in zip(noisy_pts, true_pts, quality):
        ei = ignore_quality.update(p)
        eu = use_quality.update(p, measurement_quality=q)
        e_ignore.append(_err(ei, tp))
        e_use.append(_err(eu, tp))

    burst_window = range(25, 40)
    mean_ignore = sum(e_ignore[i] for i in burst_window) / len(burst_window)
    mean_use = sum(e_use[i] for i in burst_window) / len(burst_window)
    assert mean_use < mean_ignore * 0.85


def test_default_quality_matches_static_r_behaviour():
    """measurement_quality=1.0 (domyslne) NIE powinno zmieniac R -- brak
    zewnetrznej informacji o jakosci = brak adaptacji R, tylko Q."""
    adaptive = AdaptiveKalmanFilter3D(process_var=0.01, measurement_var=5.0)
    import numpy as np
    adaptive.update((1.0, 1.0, 1.0))
    assert np.allclose(adaptive.R, adaptive._base_R)


def test_rejects_invalid_quality():
    adaptive = AdaptiveKalmanFilter3D()
    try:
        adaptive.update((0.0, 0.0, 0.0), measurement_quality=0.0)
        assert False, "powinno rzucic ValueError"
    except ValueError:
        pass
    try:
        adaptive.update((0.0, 0.0, 0.0), measurement_quality=1.5)
        assert False, "powinno rzucic ValueError"
    except ValueError:
        pass


def test_rejects_invalid_constructor_args():
    try:
        AdaptiveKalmanFilter3D(innovation_window=0)
        assert False
    except ValueError:
        pass
    try:
        AdaptiveKalmanFilter3D(q_boost_factor=0.5)
        assert False
    except ValueError:
        pass
