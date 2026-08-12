"""
visualizer.py
---------------
Wizualizacja pipeline'u FLIGHT-TRACKING-TIMDR na danych syntetycznych:
trajektoria 3D (surowy pomiar vs wygladzona przez Kalmana), krzywizna,
torsja, helikalnosc, bramkowanie (gated) i klasyfikacja manewru w
czasie.

Uzycie:
    python3 tools/visualizer.py
    -> zapisuje PNG dla kazdego oznaczonego segmentu syntetycznego do
       tools/output/<segment>.png

Wymaga matplotlib (backend Agg -- brak wyswietlacza w tym srodowisku,
zapisuje do pliku zamiast pokazywac okno).
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (rejestruje projekcje '3d')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.flight_tracker_3d import FlightTracker3D
from core.maneuver_classifier import ManeuverClassifier, ManeuverType

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")

MANEUVER_COLORS = {
    ManeuverType.STRAIGHT: "tab:blue",
    ManeuverType.TURN: "tab:orange",
    ManeuverType.SPIRAL: "tab:green",
    ManeuverType.CORKSCREW: "tab:red",
    ManeuverType.SHARP_MANEUVER: "tab:purple",
    ManeuverType.UNKNOWN: "lightgray",
}


def _load_segments(path):
    segs = {}
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        segs.setdefault(row["segment"], {"raw": [], "true": []})
        segs[row["segment"]]["raw"].append((float(row["x"]), float(row["y"]), float(row["z"])))
        segs[row["segment"]]["true"].append((float(row["x_true"]), float(row["y_true"]), float(row["z_true"])))
    return segs


def _run_pipeline(raw_points):
    tracker = FlightTracker3D(min_curvature_speed=0.3, min_curvature=1e-4)
    clf = ManeuverClassifier()
    est_points, kappas, taus, helicals, gated_flags, labels = [], [], [], [], [], []
    for x, y, z in raw_points:
        est, curv = tracker.update(x, y, z)
        label = clf.classify(curv)
        est_points.append(est)
        kappas.append(curv.curvature)
        taus.append(curv.torsion)
        helicals.append((curv.curvature ** 2 + curv.torsion ** 2) ** 0.5)
        gated_flags.append(curv.gated)
        labels.append(label)
    return est_points, kappas, taus, helicals, gated_flags, labels


def plot_segment(name, raw_points, true_points, output_path):
    est_points, kappas, taus, helicals, gated_flags, labels = _run_pipeline(raw_points)

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"FLIGHT-TRACKING-TIMDR -- segment: {name} (dane SYNTETYCZNE)", fontsize=13)

    # 1. trajektoria 3D
    ax3d = fig.add_subplot(2, 2, 1, projection="3d")
    rx, ry, rz = zip(*raw_points)
    ex, ey, ez = zip(*est_points)
    tx, ty, tz = zip(*true_points)
    ax3d.plot(tx, ty, tz, color="black", linewidth=1, linestyle="--", label="prawdziwa (referencja)")
    ax3d.scatter(rx, ry, rz, color="lightcoral", s=6, alpha=0.4, label="pomiar surowy (zaszumiony)")
    ax3d.plot(ex, ey, ez, color="tab:blue", linewidth=1.5, label="estymata Kalmana")
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_title("Trajektoria 3D")
    ax3d.legend(fontsize=7, loc="upper left")

    steps = list(range(len(raw_points)))

    # 2. krzywizna + torsja
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(steps, kappas, color="tab:blue", label="krzywizna (kappa)")
    ax2.plot(steps, taus, color="tab:orange", label="torsja (tau)", alpha=0.8)
    ax2.plot(steps, helicals, color="tab:green", label="helikalnosc (H)", alpha=0.6, linestyle=":")
    for i, g in enumerate(gated_flags):
        if g:
            ax2.axvspan(i - 0.5, i + 0.5, color="gray", alpha=0.15)
    ax2.set_xlabel("krok")
    ax2.set_title("Krzywizna / torsja / helikalnosc (szare pasy = bramkowane)")
    ax2.legend(fontsize=7)

    # 3. klasyfikacja manewru w czasie (pasek kolorow)
    ax3 = fig.add_subplot(2, 2, 3)
    for i, label in enumerate(labels):
        ax3.axvspan(i, i + 1, color=MANEUVER_COLORS[label], linewidth=0)
    ax3.set_yticks([])
    ax3.set_xlabel("krok")
    ax3.set_title("Klasyfikacja manewru w czasie")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in MANEUVER_COLORS.values()]
    ax3.legend(handles, [t.value for t in MANEUVER_COLORS.keys()], fontsize=6, loc="upper right", ncol=2)

    # 4. blad sledzenia (estymata vs prawda)
    ax4 = fig.add_subplot(2, 2, 4)
    err = [
        ((e[0] - t[0]) ** 2 + (e[1] - t[1]) ** 2 + (e[2] - t[2]) ** 2) ** 0.5
        for e, t in zip(est_points, true_points)
    ]
    raw_err = [
        ((r[0] - t[0]) ** 2 + (r[1] - t[1]) ** 2 + (r[2] - t[2]) ** 2) ** 0.5
        for r, t in zip(raw_points, true_points)
    ]
    ax4.plot(steps, raw_err, color="lightcoral", alpha=0.6, label="blad pomiaru surowego")
    ax4.plot(steps, err, color="tab:blue", label="blad estymaty Kalmana")
    ax4.set_xlabel("krok")
    ax4.set_title("Blad vs prawdziwa pozycja")
    ax4.legend(fontsize=7)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=110)
    plt.close(fig)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    segs = {}
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_flights.csv")))
    segs.update(_load_segments(os.path.join(DATA_DIR, "synthetic_maneuvers.csv")))

    for name, data in segs.items():
        out_path = os.path.join(OUTPUT_DIR, f"{name}.png")
        plot_segment(name, data["raw"], data["true"], out_path)
        print(f"zapisano {out_path}")


if __name__ == "__main__":
    main()
