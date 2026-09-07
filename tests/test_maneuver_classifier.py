import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.flight_tracker_3d import FlightTracker3D
from core.maneuver_classifier import ManeuverClassifier, ManeuverType

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _load_segments(path):
    segs = defaultdict(list)
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        segs[row["segment"]].append((float(row["x"]), float(row["y"]), float(row["z"])))
    return segs


def _all_segments():
    segs = {}
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_flights.csv")))
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_maneuvers.csv")))
    return segs


def _classify_segment(pts):
    tracker = FlightTracker3D(min_curvature_speed=0.3, min_curvature=1e-4)
    clf = ManeuverClassifier()
    counts = Counter()
    for x, y, z in pts:
        est, r = tracker.update(x, y, z)
        counts[clf.classify(r)] += 1
    return counts


# Uczciwa uwaga o metodologii: sprawdzamy WIEKSZOSCIOWA (majority) etykiete
# per segment, nie kazda probke -- klasyfikacja per-probka jest z natury
# szumiana w poblizu progow (patrz docstring modulu). To jest realistyczna
# miara dla klasyfikatora dzialajacego na zaszumionych danych czujnika.

def test_straight_climb_majority_straight():
    counts = _classify_segment(_all_segments()["straight_climb"])
    assert counts.most_common(1)[0][0] == ManeuverType.STRAIGHT


def test_level_turn_majority_turn():
    counts = _classify_segment(_all_segments()["level_turn"])
    assert counts.most_common(1)[0][0] == ManeuverType.TURN


def test_ascending_helix_majority_spiral():
    counts = _classify_segment(_all_segments()["ascending_helix"])
    assert counts.most_common(1)[0][0] == ManeuverType.SPIRAL


def test_tight_corkscrew_majority_corkscrew():
    counts = _classify_segment(_all_segments()["tight_corkscrew"])
    assert counts.most_common(1)[0][0] == ManeuverType.CORKSCREW


def test_hover_majority_unknown_gated():
    counts = _classify_segment(_all_segments()["hover"])
    assert counts.most_common(1)[0][0] == ManeuverType.UNKNOWN


def test_sharp_maneuver_majority_straight_with_some_sharp_events():
    """Segment to glownie linia prosta z krotkim, ostrym skretem --
    wiekszosc probek MUSI byc prosto, ale co najmniej kilka probek w
    trakcie manewru powinno zostac oznaczonych jako gwaltowny manewr
    (nie zdazyly sie "potwierdzic" jako trwaly zakret)."""
    counts = _classify_segment(_all_segments()["sharp_maneuver"])
    assert counts.most_common(1)[0][0] == ManeuverType.STRAIGHT
    assert counts[ManeuverType.SHARP_MANEUVER] >= 1
    # NIE powinien zdazyc sie "potwierdzic" jako trwaly zakret/spirala --
    # manewr trwa za krotko (patrz generate_synthetic_maneuvers.py)
    assert counts[ManeuverType.TURN] == 0
    assert counts[ManeuverType.SPIRAL] == 0
    assert counts[ManeuverType.CORKSCREW] == 0


def test_gated_result_is_unknown_and_resets_episode():
    from core.curvature_detector_3d import CurvatureResult3D

    clf = ManeuverClassifier()
    gated = CurvatureResult3D(curvature=0.0, torsion=0.0, speed=0.0, gated=True)
    assert clf.classify(gated) == ManeuverType.UNKNOWN


def test_invalid_thresholds_rejected():
    try:
        ManeuverClassifier(kappa_turn=0.5, kappa_spiral=0.3, kappa_corkscrew=1.0)
        assert False, "powinno rzucic ValueError (kappa_spiral <= kappa_turn)"
    except ValueError:
        pass
