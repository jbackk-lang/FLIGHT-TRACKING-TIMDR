"""
test_twist_3d_masking_and_window_limits.py -- testy dla dwóch ustaleń
znalezionych przy audycie ekosystemu TIMDR pod kątem progów liczonych z
tej samej (potencjalnie małej) próbki, którą się testuje (Pattern B):

1. NAPRAWIONE: twist_3d(robust_loo=True) (domyślne) liczy próg dla
   każdego punktu z pominięciem tego punktu (leave-one-out) - usuwa
   teoretyczne ryzyko maskowania (duży outlier zawyżający własny próg).
   Wsteczna kompatybilność: robust_loo=False odtwarza stare zachowanie
   dokładnie.

2. UCZCIWE OGRANICZENIE (ważniejsze praktycznie niż maskowanie dla tego
   modułu): na typowej długości toru (~170 punktów) z manewrem ~7%
   długości, wykrywalność jest już solidna w szerokim zakresie amplitud
   NIEZALEŻNIE od robust_loo - maskowanie nie jest tu dominującym
   problemem. Dominujący problem to CAŁKOWITA nieczułość na BARDZO
   KRÓTKICH torach (~15-25 punktów) z KRÓTKIMI manewrami (1-3 punkty),
   niezależnie od amplitudy (sprawdzone do +3200m) - okno wygładzania
   wymagane do zwalczenia wzmacniania szumu (Pattern A) jest szersze niż
   sam manewr, więc go wygładza razem z sąsiednim spokojnym lotem.
"""
import numpy as np
import pytest
from timdr_flight import TIMDRFlight


@pytest.fixture
def timdr():
    return TIMDRFlight()


def _build_maneuver_track(timdr, amp, n_pre, n_maneuver, n_post):
    R, omega, c = 4000.0, 0.04, 5.0
    n = n_pre + n_maneuver + n_post
    t = np.linspace(0, (n - 1) * 1.5, n)
    x = np.zeros(n); y = np.zeros(n); z = np.zeros(n)
    x[:n_pre] = R * np.cos(omega * t[:n_pre])
    y[:n_pre] = R * np.sin(omega * t[:n_pre])
    z[:n_pre] = 8000 + c * t[:n_pre]

    seg = slice(n_pre, n_pre + n_maneuver)
    t_seg = t[seg]
    x[seg] = R * np.cos(omega * t_seg)
    y[seg] = R * np.sin(omega * t_seg)
    period = 1.5 * max(n_maneuver - 1, 1)
    z[seg] = 8000 + c * t_seg + amp * np.sin(2 * np.pi * (t_seg - t[n_pre]) / period)

    tail = slice(n_pre + n_maneuver, n)
    x[tail] = R * np.cos(omega * t[tail])
    y[tail] = R * np.sin(omega * t[tail])
    z[tail] = z[n_pre + n_maneuver - 1] + c * (t[tail] - t[n_pre + n_maneuver - 1])

    lat, lon = timdr._inverse_local_xy(np.column_stack([x, y]), 50.0, 50.0, 20.0)
    return np.column_stack([lat, lon, z, t]), seg


class TestLeaveOneOutMasking:
    def test_robust_loo_matches_old_behavior_on_straight_flight(self, timdr):
        """Na locie prostoliniowym (brak anomalii) robust_loo=True i
        robust_loo=False powinny sie zgadzac (oba: brak falszywych
        alarmow) - naprawa nie powinna wprowadzac nowych alarmow tam,
        gdzie ich nie bylo."""
        n = 120
        t = np.linspace(0, 600, n)
        lat = np.full(n, 50.0) + 0.0005 * t / 600.0
        lon = np.full(n, 20.0) + 0.0005 * t / 600.0
        alt = np.full(n, 10000.0)
        track = np.column_stack([lat, lon, alt, t])

        idx_loo = timdr.twist_3d(track, tau_factor=3.0, robust_loo=True)
        idx_old = timdr.twist_3d(track, tau_factor=3.0, robust_loo=False)
        assert len(idx_loo) == 0
        assert len(idx_old) == 0

    def test_robust_loo_still_detects_real_maneuver(self, timdr):
        """Naprawa nie powinna pogorszyc wykrywalnosci genuine manewru
        na torze typowej dlugosci (regresja wzgledem
        test_twist_3d_wykrywa_ostry_manewr w test_frenet_serret.py)."""
        track, seg = _build_maneuver_track(timdr, amp=400, n_pre=80, n_maneuver=12, n_post=80)
        idx = timdr.twist_3d(track, tau_factor=3.0, window_seconds=6.0, robust_loo=True)
        in_seg = sum(1 for i in idx if seg.start <= i < seg.stop)
        assert in_seg > 0, "naprawa (leave-one-out) nie powinna zgubic wykrywalnosci realnego manewru"

    def test_detection_robust_across_amplitude_range_on_typical_length_track(self, timdr):
        """Na torze typowej dlugosci (~170 punktow), wykrywalnosc
        pozostaje stabilna w SZEROKIM zakresie amplitud - dokumentuje, ze
        maskowanie (Pattern B) NIE jest tu dominujacym problemem przy
        typowym rozmiarze probki (w przeciwienstwie do znalezionego
        wczesniej przypadku N~25 w GIA-TIMDR)."""
        for amp in (50, 200, 800, 1600):
            track, seg = _build_maneuver_track(timdr, amp=amp, n_pre=80, n_maneuver=12, n_post=80)
            idx = timdr.twist_3d(track, tau_factor=3.0, window_seconds=6.0, robust_loo=True)
            in_seg = sum(1 for i in idx if seg.start <= i < seg.stop)
            assert in_seg >= 10, f"amp={amp}: oczekiwano solidnej detekcji (>=10/12), uzyskano {in_seg}/12"


class TestWindowVsManeuverDurationLimit:
    """Dokumentuje (jako test, nie tylko proze) prawdziwe, dominujace
    ograniczenie znalezione przy tej naprawie: bardzo krotkie manewry na
    bardzo krotkich torach nie sa wykrywane NIEZALEZNIE OD AMPLITUDY,
    bo okno wygladzania (potrzebne przeciw Pattern A) usrednia manewr z
    sasiednim spokojnym lotem."""

    def test_short_maneuver_on_short_track_undetected_regardless_of_amplitude(self, timdr):
        detections = []
        for amp in (50, 400, 3200):
            track, seg = _build_maneuver_track(timdr, amp=amp, n_pre=10, n_maneuver=3, n_post=10)
            idx = timdr.twist_3d(track, tau_factor=3.0, window_seconds=6.0, robust_loo=True)
            in_seg = sum(1 for i in idx if seg.start <= i < seg.stop)
            detections.append(in_seg)
        assert all(d == 0 for d in detections), (
            f"oczekiwano CALKOWITEGO braku detekcji niezaleznie od amplitudy na krotkim torze "
            f"z krotkim manewrem (znane, udokumentowane ograniczenie): {detections}"
        )

    def test_same_maneuver_detected_when_track_and_maneuver_are_long_enough(self, timdr):
        """Kontrast: TEN SAM ksztalt manewru (proporcjonalnie), na
        dluzszym torze z dluzszym manewrem, JEST wykrywany - potwierdza,
        ze problem jest w relacji okno/dlugosc manewru, nie w formule."""
        track, seg = _build_maneuver_track(timdr, amp=400, n_pre=80, n_maneuver=12, n_post=80)
        idx = timdr.twist_3d(track, tau_factor=3.0, window_seconds=6.0, robust_loo=True)
        in_seg = sum(1 for i in idx if seg.start <= i < seg.stop)
        assert in_seg > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
