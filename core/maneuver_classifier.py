"""
maneuver_classifier.py
------------------------
Klasyfikuje typ manewru na podstawie krzywizny/torsji z
CurvatureDetector3D: prosto / zakret / spiralizacja / korkociag /
gwaltowny manewr / nieznany (gated).

Metodologia kalibracji (uczciwie, jak reszta repo)
----------------------------------------------------
Pojedyncza probka kappa z CurvatureDetector3D jest zbyt zaszumiona do
progowania wprost (szum GPS/barometru + wygladzanie Kalmana daja spory
rozrzut nawet dla idealnie prostego lotu). Klasyfikator dlatego liczy
medianę z ostatnich `smooth_window` (domyslnie 5) NIEBRAMKOWANYCH probek
kappa/tau i progowanie robi na tej wygladzonej wartosci -- to jest
DOKLADNIE ta sama metryka, ktora byla uzyta do kalibracji progow na
syntetycznych danych (data/synthetic_flights.csv +
data/synthetic_maneuvers.csv, 6 oznaczonych segmentow, przepuszczonych
przez FlightTracker3D). Zmierzone mediany wygladzonej kappa per segment:

    straight_climb   ~0.037   (sam szum, prawdziwe kappa=0)
    level_turn       ~0.083   (prawdziwe kappa=1/r=1/30=0.033, ale szum
                                pomiaru + wygladzanie Kalmana podnosza to)
    ascending_helix  ~0.169
    tight_corkscrew  ~0.577

Progi (kappa_turn=0.06, kappa_spiral=0.12, kappa_corkscrew=0.35) leza
miedzy tymi klastrami. TO SA WARTOSCI SKALIBROWANE NA SYNTETYCZNYM SZUMIE
GPS/barometru (sigma_xy=0.3m, sigma_z=0.5m) I NA TEJ KONKRETNEJ
GEOMETRII PROMIENI -- nie sa uniwersalna stala fizyczna. Przy innym
czujniku / innej skali manewrow wymagaja ponownej kalibracji na
realnych danych.

"Gwaltowny manewr" vs "zakret"/"spiralizacja": rozroznienie NIE jest
oparte na wielkosci kappa, tylko na CZASIE TRWANIA podwyzszonej,
wygladzonej krzywizny. Klasyfikator liczy dlugosc biezacego "epizodu"
(kolejnych probek z wygladzona kappa >= kappa_turn); dopoki epizod trwa
krocej niz `confirm_samples` (domyslnie 8), jest oznaczany jako
SHARP_MANEUVER -- klasyfikator online nie moze wiedziec z gory, czy
podwyzszona krzywizna jest poczatkiem trwalego zakretu czy tylko
krotkim manewrem. Dopiero po `confirm_samples` kolejnych probkach
epizod jest "potwierdzony" i klasyfikowany jako ZAKRET / SPIRALIZACJA /
KORKOCIAG wg wielkosci. To jest swiadomy kompromis (opoznienie detekcji
trwalego manewru o rzad `confirm_samples` probek), nie przeoczenie.

Walidacja: patrz tests/test_maneuver_classifier.py -- klasyfikator
poprawnie identyfikuje wiekszosciowa (majority) etykiete dla kazdego z
6 oznaczonych segmentow syntetycznych (dokladne liczby w teście, bo
klasyfikacja per-probka jest z natury szumiana na granicach progow).
"""

from collections import deque
from enum import Enum
from statistics import median
from typing import Deque

from core.curvature_detector_3d import CurvatureResult3D


class ManeuverType(Enum):
    STRAIGHT = "prosto"
    TURN = "zakret"
    SPIRAL = "spiralizacja"
    CORKSCREW = "korkociag"
    SHARP_MANEUVER = "gwaltowny_manewr"
    UNKNOWN = "nieznany"  # gated lub za mało historii do wygladzenia


DEFAULT_KAPPA_TURN = 0.06
DEFAULT_KAPPA_SPIRAL = 0.12
DEFAULT_KAPPA_CORKSCREW = 0.35
DEFAULT_TAU_SPIRAL = 0.4
DEFAULT_SMOOTH_WINDOW = 5
DEFAULT_CONFIRM_SAMPLES = 8


class ManeuverClassifier:
    def __init__(
        self,
        kappa_turn: float = DEFAULT_KAPPA_TURN,
        kappa_spiral: float = DEFAULT_KAPPA_SPIRAL,
        kappa_corkscrew: float = DEFAULT_KAPPA_CORKSCREW,
        tau_spiral: float = DEFAULT_TAU_SPIRAL,
        smooth_window: int = DEFAULT_SMOOTH_WINDOW,
        confirm_samples: int = DEFAULT_CONFIRM_SAMPLES,
    ):
        if kappa_turn <= 0 or kappa_spiral <= kappa_turn or kappa_corkscrew <= kappa_spiral:
            raise ValueError("wymagane: 0 < kappa_turn < kappa_spiral < kappa_corkscrew")
        self.kappa_turn = kappa_turn
        self.kappa_spiral = kappa_spiral
        self.kappa_corkscrew = kappa_corkscrew
        self.tau_spiral = tau_spiral
        self.smooth_window = smooth_window
        self.confirm_samples = confirm_samples
        self._kappa_hist: Deque[float] = deque(maxlen=smooth_window)
        self._tau_hist: Deque[float] = deque(maxlen=smooth_window)
        self._episode_len = 0
        self._last_type: "ManeuverType | None" = None

    def reset(self) -> None:
        self._kappa_hist.clear()
        self._tau_hist.clear()
        self._episode_len = 0
        self._last_type = None

    def classify(self, result: CurvatureResult3D) -> ManeuverType:
        if result.gated:
            self.reset()
            return ManeuverType.UNKNOWN

        self._kappa_hist.append(result.curvature)
        self._tau_hist.append(abs(result.torsion))

        if len(self._kappa_hist) < self.smooth_window:
            return ManeuverType.UNKNOWN  # jeszcze za mało historii do wygladzenia

        kappa = median(self._kappa_hist)
        tau = median(self._tau_hist)

        if kappa < self.kappa_turn:
            # licznik "przeciekajacy": pojedyncza probka szumu ponizej
            # progu nie zeruje od razu calego epizodu (bylby zbyt
            # kruchy przy realistycznym szumie pomiaru -- patrz
            # kalibracja w docstringu modulu), tylko go zmniejsza
            self._episode_len = max(0, self._episode_len - 2)
            if self._episode_len == 0:
                return ManeuverType.STRAIGHT
            return ManeuverType.SHARP_MANEUVER if self._episode_len < self.confirm_samples else self._last_confirmed_type()

        self._episode_len += 1

        if self._episode_len < self.confirm_samples:
            return ManeuverType.SHARP_MANEUVER

        result_type = self._classify_confirmed(kappa, tau)
        self._last_type = result_type
        return result_type

    def _classify_confirmed(self, kappa: float, tau: float) -> ManeuverType:
        if kappa >= self.kappa_corkscrew:
            return ManeuverType.CORKSCREW
        if kappa >= self.kappa_spiral:
            return ManeuverType.SPIRAL
        if tau >= self.tau_spiral:
            # kappa w paśmie "zakrętu", ale wyraźnie nie-płaski tor --
            # traktuj jako spiralę mimo niewielkiej krzywizny
            return ManeuverType.SPIRAL
        return ManeuverType.TURN

    def _last_confirmed_type(self) -> ManeuverType:
        return self._last_type if self._last_type is not None else ManeuverType.SHARP_MANEUVER
