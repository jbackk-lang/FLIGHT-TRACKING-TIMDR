"""
test_timdr_flight_trigger.py — testy timdr_flight_trigger.py.

Ten plik NIE re-weryfikuje geodezji/matematyki TIMDRFlight.twist()/
twist_3d() (już przetestowane w test_timdr_flight.py/test_frenet_serret.py)
- to nie jest robota dispatchera. Dwa rodzaje testów:

1. test_direction_twist_na_realnym_ostrym_zwrocie - JEDEN test
   integracyjny na prawdziwym TIMDRFlight (bez mockowania), używający
   TEGO SAMEGO toru co test_timdr_flight.py::test_ostry_zwrot_wykryty
   (lot na wschód, potem ostry zwrot 90st na południe) - dowód, że
   wpięcie faktycznie działa end-to-end. Lokalizacja=1 wyprowadzona
   ręcznie (projekcja equirectangular z lat0=49.7667, potem
   gradient/arctan2/unwrap/gradient - patrz analiza w PR/rozmowie).
2. Reszta testów wstrzykuje fałszywy `flight` (stub zwracający ustalony
   wynik z .twist()/.twist_3d(), ta sama struktura co
   TIMDRFlight.twist()/twist_3d()) - testujemy WYŁĄCZNIE logikę
   priorytetów/mapowania dispatchera.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from timdr_flight_trigger import TIMDRFlightTrigger, FlightTriggerType


# ----------------------------------------------------------------------
# 1) Test integracyjny na realnym TIMDRFlight
# ----------------------------------------------------------------------

def test_direction_twist_na_realnym_ostrym_zwrocie():
    """
    Ten sam tor co test_timdr_flight.py::test_ostry_zwrot_wykryty: lot na
    wschód (idx0-2), potem ostry zwrot 90st na południe (idx3-5), wysokość
    stała (1000m -> brak altitude_twist -> brak MANEUVER, sama
    DIRECTION_TWIST). Ręcznie wyprowadzone dtheta (gradient po INDEKSIE,
    nie czasie - patrz kod twist()): [0, -0.4998, -0.7854, -0.2856, 0, 0]
    rad -> |dtheta|>0.35 w indeksach 1 i 2 -> direction_twist=[1,2] ->
    lokalizacja (min) = 1.
    """
    track = [
        [50.0, 19.0, 1000, 0],
        [50.0, 19.2, 1000, 10],
        [50.0, 19.4, 1000, 20],
        [49.8, 19.4, 1000, 30],
        [49.6, 19.4, 1000, 40],
        [49.4, 19.4, 1000, 50],
    ]
    trigger = TIMDRFlightTrigger()
    result = trigger.analyze(track)

    assert result.triggered is True
    assert result.trigger_type == FlightTriggerType.DIRECTION_TWIST
    assert result.location == 1


# ----------------------------------------------------------------------
# 2) Testy priorytetów/mapowania z wstrzykniętym flight (stub)
# ----------------------------------------------------------------------

class _FakeFlight:
    """Stub o tym samym kontrakcie co TIMDRFlight: .twist() zwraca
    ustalony słownik {direction_twist, altitude_twist}, .twist_3d()
    zwraca ustaloną listę indeksów - niezależnie od danych wejściowych."""

    def __init__(self, twist_result, twist_3d_result=None):
        self._twist_result = twist_result
        self._twist_3d_result = twist_3d_result if twist_3d_result is not None else []

    def twist(self, *args, **kwargs):
        return self._twist_result

    def twist_3d(self, *args, **kwargs):
        return self._twist_3d_result


def test_priorytet_maneuver_gdy_kurs_i_wysokosc_sie_zgadzaja():
    fake = _FakeFlight({"direction_twist": [3, 7, 10], "altitude_twist": [3, 7, 9]})
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert result.trigger_type == FlightTriggerType.MANEUVER
    assert result.location == 3  # min ze wspolnych indeksow {3, 7}


def test_priorytet_altitude_nad_direction():
    fake = _FakeFlight({"direction_twist": [2], "altitude_twist": [9]})
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert result.trigger_type == FlightTriggerType.ALTITUDE_TWIST
    assert result.location == 9


def test_direction_gdy_tylko_direction():
    fake = _FakeFlight({"direction_twist": [4], "altitude_twist": []})
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert result.triggered is True
    assert result.trigger_type == FlightTriggerType.DIRECTION_TWIST
    assert result.location == 4


def test_torsion_gdy_reszta_pusta():
    fake = _FakeFlight({"direction_twist": [], "altitude_twist": []}, twist_3d_result=[6])
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert result.triggered is True
    assert result.trigger_type == FlightTriggerType.TORSION_ANOMALY
    assert result.location == 6


def test_none_gdy_wszystko_puste():
    fake = _FakeFlight({"direction_twist": [], "altitude_twist": []})
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert result.triggered is False
    assert result.trigger_type == FlightTriggerType.NONE
    assert result.location is None


def test_get_last_zwraca_ostatni_wynik():
    fake = _FakeFlight({"direction_twist": [4], "altitude_twist": []})
    trigger = TIMDRFlightTrigger(flight=fake)
    result = trigger.analyze([[0, 0, 0, 0]])
    assert trigger.get_last() is result
