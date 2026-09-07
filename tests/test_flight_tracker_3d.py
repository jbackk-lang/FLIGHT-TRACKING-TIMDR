import sys, os, csv, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.flight_tracker_3d import FlightTracker3D
from core.curvature_detector_3d import CurvatureDetector3D


def _load_segment(name):
    path = os.path.join(os.path.dirname(__file__), "..", "data", "synthetic_flights.csv")
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        lines = [l for l in f if not l.startswith("#")]
        for row in csv.DictReader(lines):
            if row["segment"] == name:
                rows.append(row)
    return rows


def test_synthetic_data_file_exists_and_has_all_segments():
    path = os.path.join(os.path.dirname(__file__), "..", "data", "synthetic_flights.csv")
    assert os.path.exists(path), "run data/generate_synthetic_flights.py first"
    with open(path, encoding="utf-8") as f:
        content = f.read()
    assert "SYNTETYCZNE DANE" in content, "synthetic-data disclaimer header missing"
    for seg in ["straight_climb", "level_turn", "ascending_helix", "hover"]:
        assert _load_segment(seg), f"missing segment {seg}"


def test_kalman_filtering_reduces_noise_vs_raw_measurement():
    """Sanity check that the Kalman stage is doing something useful: mean
    distance from the TRUE (noiseless) path should be lower after filtering
    than the raw noisy measurement, on the straight_climb segment."""
    rows = _load_segment("straight_climb")
    tracker = FlightTracker3D()

    raw_err, filtered_err = [], []
    for row in rows:
        x, y, z = float(row["x"]), float(row["y"]), float(row["z"])
        xt, yt, zt = float(row["x_true"]), float(row["y_true"]), float(row["z_true"])
        est, _ = tracker.update(x, y, z)

        raw_err.append(math.dist((x, y, z), (xt, yt, zt)))
        filtered_err.append(math.dist(est, (xt, yt, zt)))

    # skip the initial transient (first ~10 steps) where the filter is still converging
    assert sum(filtered_err[10:]) / len(filtered_err[10:]) < sum(raw_err[10:]) / len(raw_err[10:])


def test_ascending_helix_shows_higher_curvature_and_torsion_than_straight_climb():
    """Sanity ordering check using the noisy synthetic segments (run through
    the filter first, same as real usage): a genuinely curving/twisting
    path should score higher on both metrics than a straight line, even
    after noise and Kalman smoothing."""
    def mean_metrics(seg_name):
        rows = _load_segment(seg_name)
        tracker = FlightTracker3D(min_curvature_speed=0.5)
        kappas, taus = [], []
        for row in rows:
            x, y, z = float(row["x"]), float(row["y"]), float(row["z"])
            _, curv = tracker.update(x, y, z)
            if not curv.gated:
                kappas.append(curv.curvature)
                taus.append(abs(curv.torsion))
        return sum(kappas) / len(kappas), sum(taus) / len(taus)

    kappa_straight, tau_straight = mean_metrics("straight_climb")
    kappa_helix, tau_helix = mean_metrics("ascending_helix")

    assert kappa_helix > kappa_straight
    assert tau_helix > tau_straight


def test_hover_segment_mostly_gated_after_filtering():
    rows = _load_segment("hover")
    tracker = FlightTracker3D(min_curvature_speed=1.0)
    gated_count = 0
    total = 0
    for row in rows:
        x, y, z = float(row["x"]), float(row["y"]), float(row["z"])
        _, curv = tracker.update(x, y, z)
        total += 1
        if curv.gated:
            gated_count += 1
    assert gated_count / total > 0.8
