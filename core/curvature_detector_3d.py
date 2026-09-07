"""
CurvatureDetector3D -- krzywizna I torsja trajektorii w prawdziwym 3D,
dla modeli latających (dron/RC) z realnym sygnałem wysokości (barometr,
GPS altitude, itp.) -- nie z fabrykowanym pseudo-Z jak w odrzuconym
wcześniej wariancie "THE-GEO PRO 2D->3D".

Uczciwe zastrzeżenie o pochodzeniu formuł
-------------------------------------------
Oryginalny pseudokod z repo `jbackk-lang/THE` liczył torsję jako:

    tau = dot(cross(D_t-2, D_t-1), D_t) / G^2

gdzie D to znormalizowane wektory kierunku. WYGLĄDA rozsądnie (iloczyn
mieszany trzech kolejnych kierunków jako miara "wychodzenia z płaszczyzny"),
ale sprawdzone numerycznie na trajektorii o ZNANEJ analitycznie torsji
(idealna helisa 3D: x=r*cos(t), y=r*sin(t), z=c*t, gdzie krzywizna i torsja
mają dokładne wzory kappa=r/(r^2+c^2), tau=c/(r^2+c^2)) -- ta formuła
systematycznie ZBIEGA DO ZERA zamiast do prawdziwej wartości torsji, gdy
zagęszcza się próbkowanie. To jest realny błąd matematyczny w oryginalnym
wzorze, nie tylko problem z fabrykowanym Z jak poprzednio -- ujawnił się
dopiero przy rygorystycznym teście zbieżności z wzorem analitycznym,
którego nie było wcześniej (przy prawdziwych trasach GPS nie ma znanej
"prawdziwej" torsji do porównania, tylko przybliżona zmiana kursu).

Poprawna formuła (zaimplementowana tu) używa standardowych wzorów
różniczkowych na krzywiznę i torsję krzywej parametrycznej r(t), estymując
pochodne (prędkość v, przyspieszenie a, szarpnięcie/jerk j) skończonymi
różnicami z 4 kolejnych pozycji:

    kappa = |v x a| / |v|^3
    tau   = det(v, a, j) / |v x a|^2

Sprawdzone numerycznie na tej samej helisie: błąd < 0.01% już przy
dt=0.02 (patrz tests/test_curvature_detector_3d.py -- test na dokładność
analityczną). Płaski okrąg (c=0) poprawnie daje tau=0 (wykrywa płaskość),
linia prosta w dowolnym kierunku 3D daje kappa=tau=0.

Ten sam problem co w wersji 2D: przy małej prędkości (|v]|~0) dzielenie
przez |v|^3 i |v x a|^2 wzmacnia szum pomiaru do fałszywych pików. Stąd
próg min_speed poniżej którego wynik jest bramkowany (gated=True, wartości
0.0) zamiast liczony. W wersji 2D próg dobrano empirycznie z realnych
tras GPS (3.0 m/krok); tu, ponieważ repo buduje się na razie na danych
SYNTETYCZNYCH (patrz data/generate_synthetic_flights.py), próg jest
parametrem do skalibrowania na realnych logach lotu, gdy się pojawią --
domyślna wartość (1.0 jednostki/krok) jest rozsądnym punktem startowym,
nie zwalidowaną stałą.

Druga, osobna poprawka (znaleziona później, przy okazji sprawdzania
wariantu pseudokodu "THE_GEO_PRO_4D_Radar"): zabezpieczenie
`if cross_norm < 1e-12` przed torsją chroni tylko przed dosłownym
dzieleniem przez zero -- NIE chroni przed wzmocnieniem szumu, gdy
trajektoria jest niemal (ale nie dokładnie) prostoliniowa. Przykład:
v=(1,0,0), a=(1,1e-6,0) (typowy szum kierunku na "prostym" odcinku lotu)
daje cross_norm=1e-6 -- nie łapie się w próg 1e-12, więc tau wychodzi
rzędu 1e6 zamiast ~0. Poprawka: bramkowanie torsji na podstawie samej
krzywizny kappa (fizycznie sensowna wielkosc), nie na cross_norm wprost:
jeśli kappa < min_curvature -> tau = 0. Domyślne min_curvature=1e-4,
podobnie jak min_speed, to punkt startowy wymagający kalibracji na
realnych danych.
"""

from __future__ import annotations

import numpy as np
from collections import deque
from dataclasses import dataclass
from typing import Tuple


DEFAULT_MIN_SPEED = 1.0
DEFAULT_MIN_CURVATURE = 1e-4


@dataclass
class CurvatureResult3D:
    curvature: float
    torsion: float
    speed: float          # |v|, do diagnostyki i kalibracji progu
    gated: bool


class CurvatureDetector3D:
    """
    Śledzi ostatnie 4 pozycje (x, y, z) i liczy krzywiznę + torsję
    trajektorii metodą skończonych różnic (patrz docstring modułu).

    Użycie:
        det = CurvatureDetector3D(min_speed=1.0)
        for x, y, z in trajektoria_lotu:
            wynik = det.update(x, y, z)
            if not wynik.gated and wynik.torsion > próg:
                ...  # podejrzana zmiana płaszczyzny lotu (np. korkociąg)
    """

    def __init__(
        self,
        min_speed: float = DEFAULT_MIN_SPEED,
        min_curvature: float = DEFAULT_MIN_CURVATURE,
    ):
        if min_speed < 0:
            raise ValueError("min_speed nie może być ujemne")
        if min_curvature < 0:
            raise ValueError("min_curvature nie może być ujemne")
        self.min_speed = min_speed
        self.min_curvature = min_curvature
        self._positions: deque[Tuple[float, float, float]] = deque(maxlen=4)

    def reset(self) -> None:
        self._positions.clear()

    def update(self, x: float, y: float, z: float) -> CurvatureResult3D:
        self._positions.append((x, y, z))

        if len(self._positions) < 4:
            return CurvatureResult3D(curvature=0.0, torsion=0.0, speed=0.0, gated=True)

        p_t3, p_t2, p_t1, p_t = (np.array(p, dtype=float) for p in self._positions)

        v = p_t - p_t1
        a = p_t - 2 * p_t1 + p_t2
        j = p_t - 3 * p_t1 + 3 * p_t2 - p_t3

        speed = float(np.linalg.norm(v))

        if speed < self.min_speed:
            return CurvatureResult3D(curvature=0.0, torsion=0.0, speed=speed, gated=True)

        v_cross_a = np.cross(v, a)
        cross_norm = float(np.linalg.norm(v_cross_a))

        kappa = cross_norm / speed**3

        if kappa < self.min_curvature:
            # tor lokalnie prostoliniowy -- bramkujemy na podstawie kappa,
            # nie samego cross_norm==0/1e-12 (ktore chroni tylko przed
            # dzieleniem przez doslowne zero, nie przed wzmacnianiem
            # szumu gdy trajektoria jest niemal, ale nie dokladnie,
            # prosta -- patrz docstring modulu)
            tau = 0.0
        else:
            tau = float(np.dot(v_cross_a, j)) / cross_norm**2

        return CurvatureResult3D(curvature=kappa, torsion=tau, speed=speed, gated=False)
