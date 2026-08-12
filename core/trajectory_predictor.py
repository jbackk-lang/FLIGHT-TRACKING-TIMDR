"""
trajectory_predictor.py
-------------------------
Przewiduje przyszla pozycje na podstawie biezacej krzywizny/torsji
(zakladajac, ze biezacy manewr -- stale kappa, tau, predkosc -- trwa
dalej), zamiast zwyklej ekstrapolacji liniowej ze stala predkoscia.

Metoda: calkowanie rownan Freneta-Serreta (ramka T,N,B wzdluz krzywej)
metoda Rungego-Kutty 4. rzedu, z zalozeniem lokalnie stalej krzywizny i
torsji na horyzoncie predykcji -- to jest UPROSZCZENIE (prawdziwy lot
zmienia kappa/tau w czasie), analogiczne do zalozenia stalej predkosci
w zwyklym filtrze Kalmana. Dokladne przy krotkim horyzoncie predykcji i
dla manewrow o wolno zmieniajacej sie geometrii (constant-rate turn,
ustabilizowana spirala) -- NIE dla gwaltownych/nieprzewidywalnych
manewrow (do tego sluzy raczej maneuver_classifier.SHARP_MANEUVER jako
sygnal "nie ufaj predykcji").

Ramka Freneta z v (predkosc), a (przyspieszenie):
    T = v / |v|                    (styczna)
    B = (v x a) / |v x a|           (binormalna)
    N = B x T                        (normalna, dopelnia ramke prawoskretna)

Rownania Freneta-Serreta (parametryzacja dlugoscia luku s):
    dT/ds =  kappa * N
    dN/ds = -kappa * T + tau * B
    dB/ds = -tau * N
    dP/ds =  T

Walidacja: sprawdzone numerycznie na analitycznej helisie (r=5, c=1) --
przewidywana pozycja po delta_t=0.5/1.0/2.0s (przy 50 podkrokach RK4 na
predykcje) zgadza sie z pozycja analityczna z bledem < 1e-6 (patrz
tests/test_trajectory_predictor.py). Blad linowej ekstrapolacji (stala
predkosc, bez krzywizny) na tej samej helisie jest o rzedy wielkosci
wiekszy -- porownanie w tym samym teście.

Zdegenerowane przypadki: gdy tor lokalnie prostoliniowy (kappa~0, z
CurvatureDetector3D po bramkowaniu min_curvature), N i B nie sa
zdefiniowane z (v,a) -- predyktor wtedy spada do zwyklej ekstrapolacji
liniowej P + v*delta_t (poprawne zachowanie graniczne: kappa->0 oznacza
linie prosta).
"""

import math
from typing import List, Tuple

Point3 = Tuple[float, float, float]


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _scale(a, s):
    return tuple(x * s for x in a)


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def _unit(v):
    n = _norm(v)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    return tuple(x / n for x in v)


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _frenet_deriv(state, kappa, tau):
    P, T, N, B = state
    dP = T
    dT = _scale(N, kappa)
    dN = _add(_scale(T, -kappa), _scale(B, tau))
    dB = _scale(N, -tau)
    return (dP, dT, dN, dB)


def _rk4_step(state, kappa, tau, ds):
    def st_add(s1, s2, sc):
        return tuple(_add(x, _scale(y, sc)) for x, y in zip(s1, s2))

    k1 = _frenet_deriv(state, kappa, tau)
    k2 = _frenet_deriv(st_add(state, k1, ds / 2), kappa, tau)
    k3 = _frenet_deriv(st_add(state, k2, ds / 2), kappa, tau)
    k4 = _frenet_deriv(st_add(state, k3, ds), kappa, tau)
    return tuple(
        tuple(x + ds / 6 * (a + 2 * b + 2 * c + d) for x, a, b, c, d in zip(s0, K1, K2, K3, K4))
        for s0, K1, K2, K3, K4 in zip(state, k1, k2, k3, k4)
    )


def _renormalize(state):
    P, T, N, B = state
    T = _unit(T)
    B = _unit(_cross(T, N))
    N = _cross(B, T)
    return (P, T, N, B)


class TrajectoryPredictor:
    """
    Przewiduje przyszle pozycje na horyzoncie kilku krokow, zakladajac
    ze biezaca krzywizna/torsja/predkosc utrzymuja sie (patrz docstring
    modulu -- to jest uproszczenie, nie prognoza pelnej dynamiki lotu).
    """

    def __init__(self, substeps_per_step: int = 20):
        if substeps_per_step < 1:
            raise ValueError("substeps_per_step musi byc >= 1")
        self.substeps_per_step = substeps_per_step

    def predict(
        self,
        position: Point3,
        velocity: Point3,
        acceleration: Point3,
        curvature: float,
        torsion: float,
        dt: float,
        steps: int,
    ) -> List[Point3]:
        """Zwraca liste `steps` przyszlych pozycji, co `dt` sekund."""
        speed = _norm(velocity)
        if speed < 1e-9:
            # brak ruchu -- brak sensownego kierunku do ekstrapolacji
            return [position for _ in range(steps)]

        T0 = _unit(velocity)
        cross_va = _cross(velocity, acceleration)
        cross_norm = _norm(cross_va)

        if cross_norm < 1e-9 or curvature < 1e-9:
            # tor lokalnie prostoliniowy -- brak zdefiniowanej ramki
            # Freneta, spadamy do zwyklej ekstrapolacji liniowej
            return [
                _add(position, _scale(T0, speed * dt * k))
                for k in range(1, steps + 1)
            ]

        B0 = _unit(cross_va)
        N0 = _cross(B0, T0)

        state = (position, T0, N0, B0)
        results = []
        ds = (speed * dt) / self.substeps_per_step
        for _ in range(steps):
            for _ in range(self.substeps_per_step):
                state = _rk4_step(state, curvature, torsion, ds)
                state = _renormalize(state)
            results.append(state[0])
        return results
