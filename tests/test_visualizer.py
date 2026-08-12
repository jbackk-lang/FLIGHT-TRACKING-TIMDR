import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib
matplotlib.use("Agg")

from tools.visualizer import plot_segment, _load_segments, DATA_DIR


def test_plot_segment_produces_valid_png(tmp_path):
    segs = _load_segments(os.path.join(DATA_DIR, "synthetic_flights.csv"))
    data = segs["ascending_helix"]
    out_path = str(tmp_path / "test_output.png")
    plot_segment("ascending_helix", data["raw"], data["true"], out_path)
    assert os.path.exists(out_path)
    assert os.path.getsize(out_path) > 1000  # nie pusty/uszkodzony plik


def test_all_labeled_segments_load_correctly():
    segs = {}
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_flights.csv")))
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_maneuvers.csv")))
    expected = {"straight_climb", "level_turn", "ascending_helix", "hover", "tight_corkscrew", "sharp_maneuver"}
    assert set(segs.keys()) == expected
    for name, data in segs.items():
        assert len(data["raw"]) == len(data["true"])
        assert len(data["raw"]) > 0
