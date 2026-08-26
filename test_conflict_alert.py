"""
test_conflict_alert.py -- testy dla conflict_alert() (przewidywana
separacja dwoch torow), dodanego jako rozszerzenie o STCA-podobna funkcje
po sesji porownujacej TIMDR z realnymi systemami ATC (progi 5NM/1000ft to
prawdziwe minima ICAO en-route, nie wymyslone liczby - zob. README).

WAZNE: to prototyp offline na statycznych torach, NIE certyfikowany
system bezpieczenstwa. Testy sprawdzaja tylko, ze geometria/kinematyka
dziala poprawnie na jednoznacznych scenariuszach syntetycznych.
"""
import numpy as np
import pytest
from timdr_flight import TIMDRFlight


@pytest.fixture
def timdr():
    return TIMDRFlight()


def _prosty_tor(lat0, lon0, alt0, heading_deg, speed_mps, n=10, dt=10.0, t_start=0.0):
    """Prosty, jednostajny lot (stala predkosc, stala wysokosc) - do
    konstruowania scenariuszy zbieznych/rownoleglych torow."""
    t = t_start + np.arange(n) * dt
    dist = speed_mps * (t - t[0])
    brg = np.deg2rad(heading_deg)
    dx = dist * np.sin(brg)  # wschod
    dy = dist * np.cos(brg)  # polnoc
    lat = lat0 + np.rad2deg(dy / 6371000.0)
    lon = lon0 + np.rad2deg(dx / (6371000.0 * np.cos(np.deg2rad(lat0))))
    alt = np.full(n, alt0)
    return np.column_stack([lat, lon, alt, t])


def test_konflikt_wykryty_gdy_dwa_samoloty_zbiegaja_na_tej_samej_wysokosci(timdr):
    """Dwa samoloty lecace na tej samej wysokosci, kursami ktore prowadza
    wprost na siebie (glowa w glowe), blisko - powinno wykryc konflikt."""
    # samolot A: leci na wschod, zaczyna 20km na zachod od punktu spotkania
    track_a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)
    # samolot B: leci na zachod, zaczyna 20km na wschod, ta sama wysokosc
    track_b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200, n=5, dt=10.0)

    result = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0)
    assert result["conflict"] is True
    assert result["time_to_conflict_s"] is not None
    assert result["min_vertical_ft"] < 1.0  # ta sama wysokosc caly czas


def test_brak_konfliktu_gdy_duza_separacja_pionowa(timdr):
    """Te same tory poziomo (dokladnie zbiezne), ale 5000ft roznicy
    wysokosci - separacja pionowa (>1000ft) powinna wykluczyc alarm,
    mimo ze poziomo doszloby do zderzenia."""
    track_a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)
    track_b = _prosty_tor(50.0, 20.5, 11524.0, heading_deg=270, speed_mps=200, n=5, dt=10.0)  # +5000ft

    result = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0)
    assert result["conflict"] is False
    assert result["min_vertical_ft"] > 1000.0


def test_brak_konfliktu_gdy_lot_rownolegly_z_bezpieczna_separacja(timdr):
    """Dwa samoloty lecace tym samym kursem, rownolegle, w odleglosci
    10 NM (> prog 5NM) - nigdy sie nie zbliza, wiec brak alarmu."""
    track_a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)
    # 10 NM na polnoc = ok. 18520 m = ok. 0.1666 st szerokosci
    track_b = _prosty_tor(50.0 + 18520.0 / 111320.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)

    result = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0)
    assert result["conflict"] is False
    assert result["min_horizontal_nm"] > 5.0


def test_wykrywa_konflikt_w_konkretnym_oknie_czasowym(timdr):
    """Sprawdza, ze time_to_conflict_s jest sensowne (dodatnie, mniejsze
    niz lookahead) dla scenariusza z jednoznacznym momentem zblizenia."""
    track_a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)
    track_b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200, n=5, dt=10.0)

    result = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0)
    assert result["conflict"] is True
    assert 0.0 <= result["time_to_conflict_s"] <= 180.0


def test_progi_ftl_nm_konfigurowalne(timdr):
    """Zaostrzenie progu poziomego powinno moc wylaczyc alarm, ktory
    zachodzil przy progu domyslnym (potwierdza ze progi sa faktycznie
    uzywane, nie na sztywno)."""
    track_a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200, n=5, dt=10.0)
    track_b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200, n=5, dt=10.0)

    domyslny = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0, horizontal_nm=5.0)
    zaostrzony = timdr.conflict_alert(track_a, track_b, lookahead_seconds=180.0, horizontal_nm=0.001)
    assert domyslny["conflict"] is True
    assert zaostrzony["conflict"] is False
