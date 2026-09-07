"""
Generuje SYNTETYCZNE dane IMU (przyspieszeniomierz) do tych samych
segmentow lotu co synthetic_flights.csv / synthetic_maneuvers.csv --
potrzebne do core/imu_radar_fusion_3d.py.

Prawdziwe przyspieszenie liczone jest roznica centralna z KOLUMN
x_true/y_true/z_true (bezszumowa referencja, juz obecna w obu plikach
CSV) -- NIE z zaszumionych pomiarow pozycji. To jest uproszczenie:
prawdziwy IMU mierzy przyspieszenie bezposrednio (nie przez
rozniczkowanie pozycji), tutaj odtwarzamy je z generatora trajektorii,
bo to jedyne dostepne "zrodlo prawdy" w tym repo.

Model szumu przyspieszeniomierza (orientacyjny, typowy dla tanszego
MEMS, NIE zmierzony na prawdziwym sprzecie):
  - szum bialy: sigma_accel = 0.15 m/s^2 na probke
  - dryft (bias): random walk, sigma_bias_step = 0.01 m/s^2 na krok,
    inicjalizowany losowym offsetem na starcie kazdego segmentu
    (typowe dla MEMS: bias nie jest stala, "pełznie" w czasie)

IMU probkuje tu z TA SAMA CZESTOTLIWOSCIA co pomiar pozycji (radar) --
w realnym systemie IMU zwykle probkuje duzo szybciej (setki Hz vs
kilka-kilkanascie Hz radaru/GPS). To jest swiadome uproszczenie, zeby
nie komplikowac integracji w core/imu_radar_fusion_3d.py -- jawnie
udokumentowane, nie ukryte.

Uzycie:
    python3 data/generate_synthetic_imu.py
    -> zapisuje data/synthetic_imu.csv
"""
import csv
import os
import random
from collections import defaultdict

random.seed(44)

SIGMA_ACCEL = 0.15
SIGMA_BIAS_STEP = 0.01
SIGMA_BIAS_INIT = 0.1


def _load_true_trajectory(path):
    segs = defaultdict(list)
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        segs[row["segment"]].append(
            (float(row["x_true"]), float(row["y_true"]), float(row["z_true"]))
        )
    return segs


def _central_diff_accel(pts):
    """Przyspieszenie = druga pochodna pozycji, roznica centralna
    (dt=1 krok, spojne z reszta repo -- patrz zastrzezenie w
    kalman_filter_3d.py o F/Q budowanych dla dt=1)."""
    n = len(pts)
    accel = []
    for i in range(n):
        if i == 0 or i == n - 1:
            accel.append((0.0, 0.0, 0.0))
        else:
            ax = pts[i + 1][0] - 2 * pts[i][0] + pts[i - 1][0]
            ay = pts[i + 1][1] - 2 * pts[i][1] + pts[i - 1][1]
            az = pts[i + 1][2] - 2 * pts[i][2] + pts[i - 1][2]
            accel.append((ax, ay, az))
    return accel


def main():
    data_dir = os.path.dirname(__file__)
    segs = {}
    segs.update(_load_true_trajectory(os.path.join(data_dir, "synthetic_flights.csv")))
    segs.update(_load_true_trajectory(os.path.join(data_dir, "synthetic_maneuvers.csv")))

    rows = []
    for seg_name, pts in segs.items():
        true_accel = _central_diff_accel(pts)
        bias = [random.gauss(0, SIGMA_BIAS_INIT) for _ in range(3)]
        for i, (ax, ay, az) in enumerate(true_accel):
            bias = [b + random.gauss(0, SIGMA_BIAS_STEP) for b in bias]
            meas = (
                ax + bias[0] + random.gauss(0, SIGMA_ACCEL),
                ay + bias[1] + random.gauss(0, SIGMA_ACCEL),
                az + bias[2] + random.gauss(0, SIGMA_ACCEL),
            )
            rows.append({
                "segment": seg_name,
                "step": i,
                "ax_true": round(ax, 5),
                "ay_true": round(ay, 5),
                "az_true": round(az, 5),
                "ax_meas": round(meas[0], 5),
                "ay_meas": round(meas[1], 5),
                "az_meas": round(meas[2], 5),
            })

    out_path = os.path.join(data_dir, "synthetic_imu.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        f.write("# SYNTETYCZNE DANE IMU -- nie sa to pomiary z prawdziwego czujnika.\n")
        f.write("# Wygenerowane przez generate_synthetic_imu.py z roznicy centralnej\n")
        f.write("# trajektorii x_true/y_true/z_true + szum bialy (sigma=0.15 m/s^2) +\n")
        f.write("# dryft/bias (random walk, sigma_step=0.01 m/s^2). Patrz docstring.\n")
        writer = csv.DictWriter(f, fieldnames=["segment", "step", "ax_true", "ay_true", "az_true", "ax_meas", "ay_meas", "az_meas"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"zapisano {len(rows)} wierszy do {out_path}")


if __name__ == "__main__":
    main()
