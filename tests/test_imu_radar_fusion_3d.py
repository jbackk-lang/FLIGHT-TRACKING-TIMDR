import csv
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.kalman_filter_3d import KalmanFilter3D
from core.imu_radar_fusion_3d import ImuRadarFusion3D, ImuDeadReckoning3D

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _err(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _load(path, seg, cols):
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    return [tuple(float(r[c]) for c in cols) for r in reader if r["segment"] == seg]


def _load_segment(seg):
    radar = _load(os.path.join(DATA_DIR, "synthetic_maneuvers.csv"), seg, ["x", "y", "z"])
    true = _load(os.path.join(DATA_DIR, "synthetic_maneuvers.csv"), seg, ["x_true", "y_true", "z_true"])
    imu = _load(os.path.join(DATA_DIR, "synthetic_imu.csv"), seg, ["ax_meas", "ay_meas", "az_meas"])
    assert len(radar) == len(true) == len(imu)
    return radar, true, imu


def test_fusion_beats_radar_only_during_maneuver_with_matched_measurement_var():
    """Uczciwe porownanie: TA SAMA measurement_var dla obu wariantow
    (dopasowana do rzeczywistego szumu symulacji pozycji -- patrz
    docstring modulu), zeby izolowac wklad samego IMU, nie efekt
    przypadkowo lepiej dobranego R."""
    radar_pts, true_pts, imu_meas = _load_segment("sharp_maneuver")
    MV = 0.3

    radar_only = KalmanFilter3D(process_var=0.01, measurement_var=MV, initial_position=radar_pts[0])
    fusion = ImuRadarFusion3D(process_var=0.001, measurement_var=MV, initial_position=radar_pts[0])

    radar_errs, fusion_errs = [], []
    for rp, tp, ac in zip(radar_pts, true_pts, imu_meas):
        radar_errs.append(_err(radar_only.update(rp), tp))
        fusion_errs.append(_err(fusion.update(ac, rp), tp))

    maneuver_window = range(38, 55)
    radar_mean = sum(radar_errs[i] for i in maneuver_window) / len(maneuver_window)
    fusion_mean = sum(fusion_errs[i] for i in maneuver_window) / len(maneuver_window)
    assert fusion_mean < radar_mean * 0.85  # co najmniej 15% lepiej podczas manewru


def test_imu_only_drifts_without_correction():
    """Martwa nawigacja z samego IMU (bez radaru) MUSI dryfowac znaczaco
    -- to jest uzasadnienie, dlaczego fuzja (nie sam IMU) jest potrzebna."""
    radar_pts = _load(os.path.join(DATA_DIR, "synthetic_flights.csv"), "straight_climb", ["x", "y", "z"])
    true_pts = _load(os.path.join(DATA_DIR, "synthetic_flights.csv"), "straight_climb", ["x_true", "y_true", "z_true"])
    imu_meas = _load(os.path.join(DATA_DIR, "synthetic_imu.csv"), "straight_climb", ["ax_meas", "ay_meas", "az_meas"])

    dead_reck = ImuDeadReckoning3D(initial_position=radar_pts[0])
    errs = []
    for tp, ac in zip(true_pts, imu_meas):
        errs.append(_err(dead_reck.update(ac), tp))

    # blad musi rosnac wyraznie w czasie (dryft), nie zostawac ograniczony
    early_err = sum(errs[:5]) / 5
    late_err = sum(errs[-5:]) / 5
    assert late_err > early_err * 5


def test_fusion_tracks_far_better_than_dead_reckoning_long_term():
    radar_pts = _load(os.path.join(DATA_DIR, "synthetic_flights.csv"), "ascending_helix", ["x", "y", "z"])
    true_pts = _load(os.path.join(DATA_DIR, "synthetic_flights.csv"), "ascending_helix", ["x_true", "y_true", "z_true"])
    imu_meas = _load(os.path.join(DATA_DIR, "synthetic_imu.csv"), "ascending_helix", ["ax_meas", "ay_meas", "az_meas"])

    fusion = ImuRadarFusion3D(process_var=0.001, measurement_var=0.3, initial_position=radar_pts[0])
    dead_reck = ImuDeadReckoning3D(initial_position=radar_pts[0])

    fusion_errs, dr_errs = [], []
    for rp, tp, ac in zip(radar_pts, true_pts, imu_meas):
        fusion_errs.append(_err(fusion.update(ac, rp), tp))
        dr_errs.append(_err(dead_reck.update(ac), tp))

    assert fusion_errs[-1] < dr_errs[-1] / 10  # fuzja rzedy wielkosci lepsza dlugoterminowo


def test_rejects_invalid_dt():
    try:
        ImuRadarFusion3D(dt=0.0)
        assert False, "powinno rzucic ValueError"
    except ValueError:
        pass
    try:
        ImuRadarFusion3D(dt=-1.0)
        assert False
    except ValueError:
        pass


def test_predict_then_correct_matches_update():
    """update() musi byc rownowazne predict() nastepnie correct()."""
    f1 = ImuRadarFusion3D(initial_position=(1.0, 2.0, 3.0))
    f2 = ImuRadarFusion3D(initial_position=(1.0, 2.0, 3.0))

    accel = (0.5, -0.2, 0.1)
    radar_pos = (2.0, 2.5, 3.5)

    r1 = f1.update(accel, radar_pos)
    f2.predict(accel)
    r2 = f2.correct(radar_pos)

    assert all(abs(a - b) < 1e-9 for a, b in zip(r1, r2))
