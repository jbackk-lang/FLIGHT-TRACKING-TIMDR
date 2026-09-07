"""
Generuje SYNTETYCZNE trajektorie do kalibracji i testow
maneuver_classifier.py. Rozszerza cztery segmenty z
generate_synthetic_flights.py o dwa nowe, potrzebne do rozroznienia
wszystkich 5 kategorii manewru:

  - straight_climb   (juz w synthetic_flights.csv) -> "prosto"
  - level_turn       (juz w synthetic_flights.csv) -> "zakret"
  - ascending_helix  (juz w synthetic_flights.csv) -> "spiralizacja"
  - tight_corkscrew  (NOWY, ten plik)               -> "korkociag"
  - sharp_maneuver   (NOWY, ten plik)                -> "gwaltowny_manewr"

tight_corkscrew: helisa o duzo mniejszym promieniu i wiekszej predkosci
katowej niz ascending_helix -- fizycznie "korkociag" to ciasny, szybki
obrot, w odroznieniu od lagodnego wznoszenia spiralnego.

sharp_maneuver: lot prosto, potem nagla zmiana kierunku o ~70 stopni w
ciagu 3 krokow, potem znowu prosto -- krotki, gwaltowny SKOK krzywizny,
w odroznieniu od zakretu, ktory jest utrzymany przez wiele krokow.

Te dane sa jawnie oznaczone jako syntetyczne (naglowek CSV + kolumna
'segment'), tak jak w generate_synthetic_flights.py -- nie sa to
prawdziwe logi lotu.
"""
import csv
import math
import os
import random

random.seed(43)

SIGMA_XY = 0.3
SIGMA_Z = 0.5


def noisy(p):
    x, y, z = p
    return (
        x + random.gauss(0, SIGMA_XY),
        y + random.gauss(0, SIGMA_XY),
        z + random.gauss(0, SIGMA_Z),
    )


def gen_tight_corkscrew(n=150, r=3.0, c=1.0, dt=0.15):
    return [(r * math.cos(t * dt), r * math.sin(t * dt), c * t * dt) for t in range(n)]


def gen_sharp_maneuver(n_straight1=40, n_turn=3, n_straight2=40, speed=2.5, turn_deg=70):
    """Prosto -> nagly skret o turn_deg w ciagu n_turn krokow -> prosto."""
    pts = []
    heading = 0.0
    pos = (0.0, 0.0, 30.0)
    for _ in range(n_straight1):
        dx = speed * math.cos(heading)
        dy = speed * math.sin(heading)
        pos = (pos[0] + dx, pos[1] + dy, pos[2])
        pts.append(pos)
    turn_step = math.radians(turn_deg) / n_turn
    for _ in range(n_turn):
        heading += turn_step
        dx = speed * math.cos(heading)
        dy = speed * math.sin(heading)
        pos = (pos[0] + dx, pos[1] + dy, pos[2])
        pts.append(pos)
    for _ in range(n_straight2):
        dx = speed * math.cos(heading)
        dy = speed * math.sin(heading)
        pos = (pos[0] + dx, pos[1] + dy, pos[2])
        pts.append(pos)
    return pts


def main():
    rows = []
    segments = {
        "tight_corkscrew": gen_tight_corkscrew(),
        "sharp_maneuver": gen_sharp_maneuver(),
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

    out_path = os.path.join(os.path.dirname(__file__), "synthetic_maneuvers.csv")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        f.write("# SYNTETYCZNE DANE -- nie sa to prawdziwe logi lotu. Wygenerowane przez\n")
        f.write("# generate_synthetic_maneuvers.py do kalibracji/testow maneuver_classifier.py.\n")
        writer = csv.DictWriter(f, fieldnames=["segment", "step", "x", "y", "z", "x_true", "y_true", "z_true"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"zapisano {len(rows)} wierszy do {out_path}")


if __name__ == "__main__":
    main()
