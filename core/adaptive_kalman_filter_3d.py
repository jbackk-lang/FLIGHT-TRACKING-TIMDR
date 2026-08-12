"""
adaptive_kalman_filter_3d.py
------------------------------
Rozszerza KalmanFilter3D (stala predkosc, stale Q/R) o:

  1. Adaptacyjne Q -- proces "innovation-based adaptive estimation"
     (IAE, klasyczna technika, Mehra 1970s): sledzi energie ostatnich
     innowacji (roznica pomiar-predykcja) w oknie kroczacym. Gdy
     srednia energia innowacji przekracza `innovation_threshold` razy
     energie oczekiwana pod hipoteza "brak manewru" (slad R), Q jest
     mnozone przez `q_boost_factor` -- filtr zaczyna bardziej ufac
     nowym pomiarom niz wlasnemu modelowi stalej predkosci, bo
     systematycznie duze innowacje oznaczaja, ze cel faktycznie
     przyspiesza/skreca, nie ze to szum. Samowystarczalne -- nie
     wymaga zewnetrznego sygnalu krzywizny.

  2. Opcjonalne, zewnetrznie podawane `measurement_quality` (0,1] przy
     kazdym wywolaniu update() -- skaluje R w gore, gdy jakosc pomiaru
     jest niska (np. slabe SNR radaru). W przeciwienstwie do Q, NIE ma
     tu samowystarczalnego, wewnetrznego sposobu ocenic "jakosc
     pomiaru" bez dodatkowych danych od czujnika (np. SNR) -- domyslna
     wartosc 1.0 (brak informacji o jakosci) daje dokladnie taki sam
     wynik jak statyczne R w KalmanFilter3D. To jest swiadomie
     "martwe" bez prawdziwego zrodla danych o jakosci sygnalu -- ale
     PRZETESTOWANE z jawnie wstrzykiwana syntetyczna jakoscia (patrz
     tests/test_adaptive_kalman_filter_3d.py), zeby udowodnic, ze
     mechanizm dziala, gdy dane beda dostepne.

Walidacja: na syntetycznej trajektorii prosto->ostry skret->prosto
(taka jak data/synthetic_maneuvers.csv, segment sharp_maneuver),
adaptacyjne Q daje mniejszy blad sledzenia W TRAKCIE i BEZPOSREDNIO PO
manewrze niz statyczne Q (filtr "nadraza" szybciej), przy porownywalnym
bledzie na odcinkach prostych (patrz test). To jest typowy kompromis
adaptacyjnego Kalmana: szybsza reakcja na manewr kosztem odrobine
wiekszego szumu na spokojnych odcinkach.
"""
from collections import deque

import numpy as np

from core.kalman_filter_3d import KalmanFilter3D

DEFAULT_INNOVATION_WINDOW = 3
DEFAULT_Q_BOOST_FACTOR = 50.0
DEFAULT_INNOVATION_THRESHOLD = 1.3


class AdaptiveKalmanFilter3D(KalmanFilter3D):
    def __init__(
        self,
        process_var: float = 0.01,
        measurement_var: float = 5.0,
        initial_position=(0.0, 0.0, 0.0),
        initial_uncertainty: float = 100.0,
        innovation_window: int = DEFAULT_INNOVATION_WINDOW,
        q_boost_factor: float = DEFAULT_Q_BOOST_FACTOR,
        innovation_threshold: float = DEFAULT_INNOVATION_THRESHOLD,
    ):
        super().__init__(process_var, measurement_var, initial_position, initial_uncertainty)
        if innovation_window < 1:
            raise ValueError("innovation_window musi byc >= 1")
        if q_boost_factor < 1.0:
            raise ValueError("q_boost_factor musi byc >= 1.0 (to wzmocnienie, nie tlumienie)")
        self._base_Q = self.Q.copy()
        self._base_R = self.R.copy()
        self._innovations: deque = deque(maxlen=innovation_window)
        self.q_boost_factor = q_boost_factor
        self.innovation_threshold = innovation_threshold
        self.last_maneuvering = False  # diagnostyka: czy ostatni krok uznano za manewr

    def _adapt_Q(self) -> np.ndarray:
        if len(self._innovations) < self._innovations.maxlen:
            self.last_maneuvering = False
            return self._base_Q

        innov = np.array(self._innovations)  # (window, 3)
        mean_energy = float(np.mean(np.sum(innov ** 2, axis=1)))
        baseline_energy = float(np.trace(self._base_R))
        ratio = mean_energy / max(baseline_energy, 1e-9)

        self.last_maneuvering = ratio > self.innovation_threshold
        return self._base_Q * self.q_boost_factor if self.last_maneuvering else self._base_Q

    def update(self, z, measurement_quality: float = 1.0):
        if not (0.0 < measurement_quality <= 1.0):
            raise ValueError("measurement_quality musi byc w (0, 1]")

        zx, zy, zz = z
        z_vec = np.array([[zx], [zy], [zz]])

        self.Q = self._adapt_Q()
        R_eff = self._base_R / measurement_quality  # nizsza jakosc -> wiekszy efektywny szum pomiaru

        # predykcja
        self.x = self.F @ self.x
        self.P = self.F @ self.P @ self.F.T + self.Q

        # korekcja
        y = z_vec - self.H @ self.x
        S = self.H @ self.P @ self.H.T + R_eff
        K = self.P @ self.H.T @ np.linalg.inv(S)

        self.x = self.x + K @ y
        self.P = (np.eye(6) - K @ self.H) @ self.P

        self._innovations.append(y.flatten())
        self.R = R_eff  # dla diagnostyki/introspekcji

        return (float(self.x[0, 0]), float(self.x[1, 0]), float(self.x[2, 0]))
