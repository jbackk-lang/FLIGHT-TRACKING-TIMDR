"""
Generuje SYNTETYCZNE trajektorie lotu 3D do testów sanity CurvatureDetector3D
i FlightTracker3D. Nie są to prawdziwe logi lotu -- oznaczone jawnie w
kolumnie `segment` i w nagłówku pliku CSV.

Cztery segmenty, każdy z dodanym realistycznym szumem pomiaru GPS/barometru
(sigma_xy=0.3m, sigma_z=0.5m -- typowe rzędy wielkości dla konsumenckiego
GPS+barometru, nie zmierzone, tylko orientacyjne):

1. straight_climb  -- linia prosta ze stałym wznoszeniem (kappa=tau=0)
2. level_turn       -- płaski zakręt na stałej wysokości (tau=0, kappa=1/r)
3. ascending_helix  -- wznosząca się spirala (kappa i tau ze wzoru
                        analitycznego: kappa=r/(r^2+c^2), tau=c/(r^2+c^2))
4. hover            -- prawie-postój z samym szumem czujnika (do testu
                        bramkowania przy niskiej prędkości)

Użycie:
    python3 data/generate_synthetic_flights.py
    -> zapisuje data/synthetic_flights.csv
"""
import csv
import math
import os
import random

random.seed(42)

SIGMA_XY = 0.3
SIGMA_Z = 0.5


def noisy(p):
    x, y, z = p
    return (
        x + random.gauss(0, SIGMA_XY),
        y + random.gauss(0, SIGMA_XY),
        z + random.gauss(0, SIGMA_Z),
    )


def gen_straight_climb(n=60):
    return [(2.0 * t, 1.0 * t, 0.5 * t) for t in range(n)]


def gen_level_turn(n=120, r=30.0, dt=0.05):
    return [(r * math.cos(t * dt), r * math.sin(t * dt), 50.0) for t in range(n)]


def gen_ascending_helix(n=200, r=20.0, c=3.0, dt=0.05):
    return [(r * math.cos(t * dt), r * math.sin(t * dt), c * t * dt) for t in range(n)]


def gen_hover(n=40, z0=40.0):
    return [(0.0, 0.0, z0) for _ in range(n)]


def main():
    rows = []
    segments = {
        "straight_climb": gen_straight_climb(),
        "level_turn": gen_level_turn(),
        "ascending_helix": gen_ascending_helix(),
        "hover": gen_hover(),
    }
    for seg_name, pts in segments.items():
        for i, p in enumerate(pts):
            nx, ny, nz = noisy(p)
            rows.append({
                "segment": seg_name,
                "step": i,
                "x": round(nx, 4),
                "y": round(ny, 4),
                "z": round(nz, 4),
                "x_true": round(p[0], 4),
                "y_true": round(p[1], 4),
                "z_true": round(p[2], 4),
            })

    out_path = os.path.join(os.path.dirname(__file__), "synthetic_flights.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        f.write("# SYNTETYCZNE DANE -- nie sa to prawdziwe logi lotu. Wygenerowane przez\n")
        f.write("# generate_synthetic_flights.py do testow sanity, sigma_xy=0.3m sigma_z=0.5m.\n")
        writer = csv.DictWriter(f, fieldnames=["segment", "step", "x", "y", "z", "x_true", "y_true", "z_true"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"zapisano {len(rows)} wierszy do {out_path}")


if __name__ == "__main__":
    main()
