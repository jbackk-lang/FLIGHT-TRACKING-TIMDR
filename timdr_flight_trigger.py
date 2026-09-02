# ============================================
# TIMDR Flight Trigger Module
# ============================================
#
# ROLA: czujnik integralności JEDNEGO toru lotu — NIE model, NIE predyktor
# trajektorii (do tego służy TIMDRFlight.predict()). Dispatcher nad już
# przetestowanym TIMDRFlight.twist()/twist_3d() (test_timdr_flight.py,
# test_frenet_serret.py) — jedyna jego robota: powiedzieć, KTÓRY typ
# manewru/anomalii się odpalił i GDZIE (indeks punktu w torze).
#
# To jest odpowiednik dla JEDNEGO toru tego, co conflict_alert()/
# conflict_alert_fleet() już robi dla PARY/FLOTY torów (separacja
# horyzontalna/wertykalna) — świadomie NIE duplikuje tamtej, już w pełni
# ugruntowanej logiki, tylko uzupełnia ją o integralność pojedynczego
# toru (kurs/wysokość/torsja), czego conflict_alert() nie sprawdza.
#
# Priorytet: MANEUVER (zmiana kursu I prędkości pionowej w TYM SAMYM
# punkcie — najsilniejszy dowód gwałtownego manewru, dwa niezależne
# detektory się zgadzają) > ALTITUDE_TWIST (samo strome wznoszenie/
# zniżanie — fizycznie najbardziej dotkliwe pojedynczo: ograniczenia
# strukturalne/komfortu, względy TCAS) > DIRECTION_TWIST (sama nagła
# zmiana kursu) > TORSION_ANOMALY (twist_3d — NAJMNIEJ zweryfikowany
# sygnał w tym module, patrz zastrzeżenie w docstringu frenet_serret()/
# twist_3d(): "to NIE jest zwalidowane na prawdziwych danych ADS-B") >
# NONE. Silniejszy/łączny dowód wygrywa niezależnie od tego, co jest
# chronologicznie pierwsze w torze — ta sama zasada co w reszcie
# ekosystemu TIMDR.

from enum import Enum

from timdr_flight import TIMDRFlight


class FlightTriggerType(Enum):
    MANEUVER = "combined_maneuver"
    ALTITUDE_TWIST = "altitude_twist"
    DIRECTION_TWIST = "direction_twist"
    TORSION_ANOMALY = "torsion_anomaly"
    NONE = "none"


class FlightTriggerResult:
    def __init__(self, triggered=False, trigger_type=FlightTriggerType.NONE,
                 location=None, message=""):
        self.triggered = triggered
        self.trigger_type = trigger_type
        self.location = location
        self.message = message

    def as_dict(self):
        return {
            "triggered": self.triggered,
            "type": self.trigger_type.value,
            "location": self.location,
            "message": self.message,
        }


class TIMDRFlightTrigger:
    """
    Dispatcher nad TIMDRFlight.twist()/twist_3d(). `flight` można
    wstrzyknąć (np. w testach) - domyślnie tworzy prawdziwy TIMDRFlight().
    Progi (angle_thresh, climb_rate_thresh_mps, tau_factor) to te same
    punkty startowe do dostrojenia co w reszcie ekosystemu TIMDR, nie
    wartości uniwersalne - patrz docstringi odpowiadających metod w
    timdr_flight.py.
    """

    def __init__(self, angle_thresh=0.35, climb_rate_thresh_mps=15.0,
                 tau_factor=3.0, flight=None):
        self.flight = flight if flight is not None else TIMDRFlight()
        self.angle_thresh = angle_thresh
        self.climb_rate_thresh_mps = climb_rate_thresh_mps
        self.tau_factor = tau_factor
        self.last_result = FlightTriggerResult()

    def analyze(self, track):
        tw = self.flight.twist(
            track,
            angle_thresh=self.angle_thresh,
            climb_rate_thresh_mps=self.climb_rate_thresh_mps,
        )
        dir_idx = set(int(i) for i in tw["direction_twist"])
        alt_idx = set(int(i) for i in tw["altitude_twist"])

        both = dir_idx & alt_idx
        if both:
            loc = min(both)
            return self._set_result(
                True, FlightTriggerType.MANEUVER, loc,
                "Jednoczesna nagła zmiana kursu i prędkości pionowej."
            )

        if alt_idx:
            loc = min(alt_idx)
            return self._set_result(
                True, FlightTriggerType.ALTITUDE_TWIST, loc,
                "Nagła zmiana prędkości pionowej (strome wznoszenie/zniżanie)."
            )

        if dir_idx:
            loc = min(dir_idx)
            return self._set_result(
                True, FlightTriggerType.DIRECTION_TWIST, loc,
                "Nagła zmiana kursu."
            )

        torsion_idx = self.flight.twist_3d(track, tau_factor=self.tau_factor)
        if len(torsion_idx):
            loc = int(min(torsion_idx))
            return self._set_result(
                True, FlightTriggerType.TORSION_ANOMALY, loc,
                "Anomalia skręcenia (torsji) trajektorii - sygnał "
                "niezweryfikowany na realnych danych ADS-B."
            )

        return self._set_result(
            False, FlightTriggerType.NONE, None,
            "Brak wykrytego zdarzenia integralności toru."
        )

    def _set_result(self, triggered, trigger_type, location, message):
        self.last_result = FlightTriggerResult(triggered, trigger_type, location, message)
        return self.last_result

    def get_last(self):
        return self.last_result
