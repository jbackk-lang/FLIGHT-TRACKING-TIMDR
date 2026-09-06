"""
test_noise_amplification_fix.py
--------------------------------
Test regresyjny dla NAPRAWIONEGO problemu "wzmacnianie szumu przy
potrojnym roznicowaniu" (Pattern A, patrz komentarz w
the_geo_pro_4d.py nad THE_GEO_PRO_4D_sequence()).

Metoda identyczna z tym, jak FLIGHT-TRACKING-TIMDR zmierzylo ten sam
problem na wlasnym kodzie: zaszumiona helisa (znana analitycznie
tau), porownanie odchylenia standardowego bledu tau miedzy surowym
roznicowaniem (smooth=False) a wygladzonym (smooth=True, domyslne).
"""
import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from the_geo_pro_4d import THE_GEO_PRO_4D_sequence


def noisy_helix_points(n, r=5.0, c=1.0, dt=0.05, t0=2.0, noise_std=0.05, seed=0):
    rng = random.Random(seed)
    pts = []
    for k in range(n):
        t = t0 + k * dt
        x = r * math.cos(t) + rng.gauss(0, noise_std)
        y = r * math.sin(t) + rng.gauss(0, noise_std)
        z = c * t + rng.gauss(0, noise_std)
        pts.append((x, y, z))
    return pts


def analytical_tau(r=5.0, c=1.0):
    return c / (r * r + c * c)


class TestNoiseAmplificationFix(unittest.TestCase):
    """
    NAPRAWIONE: THE_GEO_PRO_4D_sequence(smooth=True) (domyslne)
    powinno miec ISTOTNIE mniejszy blad/rozrzut tau wzgledem wartosci
    analitycznej niz smooth=False (stare, surowe zachowanie) na
    zaszumionej helisie -- ten sam wzorzec pomiaru, ktorego uzylo
    FLIGHT-TRACKING-TIMDR (tam: std(tau)~2.6 surowe vs ~0.065
    prawdziwe, 40x wzmocnienie szumu).
    """

    def test_smoothed_sequence_has_lower_tau_error_than_raw(self):
        r, c, dt = 5.0, 1.0, 0.05
        tau_true = analytical_tau(r, c)
        pts = noisy_helix_points(n=60, r=r, c=c, dt=dt, noise_std=0.05, seed=42)

        raw_results = THE_GEO_PRO_4D_sequence(pts, min_speed=0.0, smooth=False)
        smoothed_results = THE_GEO_PRO_4D_sequence(pts, min_speed=0.0, smooth=True, window=11, polyorder=3)

        raw_taus = [r_["torsion"] for r_ in raw_results if not r_["gated"]]
        smoothed_taus = [r_["torsion"] for r_ in smoothed_results if not r_["gated"]]

        raw_errors = [abs(t - tau_true) for t in raw_taus]
        smoothed_errors = [abs(t - tau_true) for t in smoothed_taus]

        raw_mean_err = sum(raw_errors) / len(raw_errors)
        smoothed_mean_err = sum(smoothed_errors) / len(smoothed_errors)

        self.assertLess(
            smoothed_mean_err, raw_mean_err * 0.7,
            f"oczekiwano wyraznej redukcji bledu tau po wygladzeniu: "
            f"surowe={raw_mean_err:.4f}, wygladzone={smoothed_mean_err:.4f} (prawdziwe tau={tau_true:.4f})"
        )

    def test_smoothing_does_not_break_clean_analytical_helix(self):
        """Sanity: na CZYSTEJ (bezszumowej) helisie wygladzanie NIE
        powinno pogorszyc zbieznosci istotnie ponizej tego, co juz
        dawalo surowe roznicowanie (TestHelixConvergence w
        test_the_geo_pro_4d.py) -- wygladzanie ma pomagac na szumie,
        nie psuc dokladnosc tam, gdzie szumu nie ma."""
        r, c, dt = 5.0, 1.0, 0.02
        tau_true = analytical_tau(r, c)
        pts = noisy_helix_points(n=30, r=r, c=c, dt=dt, noise_std=0.0, seed=1)
        results = THE_GEO_PRO_4D_sequence(pts, min_speed=0.0, smooth=True, window=9, polyorder=3)
        taus = [res["torsion"] for res in results if not res["gated"]]
        mean_err = sum(abs(t - tau_true) for t in taus) / len(taus)
        self.assertLess(mean_err / tau_true, 0.05, "blad wzgledny powinien pozostac maly na czystych danych")

    def test_smooth_false_reproduces_old_raw_behavior(self):
        """smooth=False powinno dawac DOKLADNIE to samo, co reczne
        wywolanie THE_GEO_PRO_4D na kolejnych oknach 4 punktow (stare
        zachowanie, wsteczna kompatybilnosc)."""
        from the_geo_pro_4d import THE_GEO_PRO_4D
        pts = noisy_helix_points(n=10, noise_std=0.02, seed=7)
        seq_results = THE_GEO_PRO_4D_sequence(pts, min_speed=0.0, smooth=False)
        manual_results = [
            THE_GEO_PRO_4D(pts[i], pts[i-1], pts[i-2], pts[i-3], min_speed=0.0)
            for i in range(3, len(pts))
        ]
        for a, b in zip(seq_results, manual_results):
            self.assertEqual(a, b)

    def test_polyorder_below_3_raises_clear_error_not_indexerror(self):
        pts = noisy_helix_points(n=10, noise_std=0.0, seed=1)
        with self.assertRaises(ValueError):
            THE_GEO_PRO_4D_sequence(pts, smooth=True, polyorder=2)


if __name__ == "__main__":
    unittest.main()
