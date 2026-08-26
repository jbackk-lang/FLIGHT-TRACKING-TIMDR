"""
test_conflict_fleet_airspace.py -- testy dla AIRSPACE_PRESETS i
conflict_alert_fleet(), dodanych po rekomendacji z sesji porownawczej z
realnym STCA: presety progow per typ przestrzeni (ICAO Doc 4444/ECAC -
en-route 5NM, TMA 3NM, final approach 2.5NM, wszystkie 1000ft pionowo) i
skan wszystkich par torow naraz zamiast tylko jednej pary.
"""
import numpy as np
import pytest
from timdr_flight import TIMDRFlight


@pytest.fixture
def timdr():
    return TIMDRFlight()


def _prosty_tor(lat0, lon0, alt0, heading_deg, speed_mps, n=5, dt=10.0):
    t = np.arange(n) * dt
    dist = speed_mps * t
    brg = np.deg2rad(heading_deg)
    dx = dist * np.sin(brg)
    dy = dist * np.cos(brg)
    lat = lat0 + np.rad2deg(dy / 6371000.0)
    lon = lon0 + np.rad2deg(dx / (6371000.0 * np.cos(np.deg2rad(lat0))))
    alt = np.full(n, alt0)
    return np.column_stack([lat, lon, alt, t])


def _rownolegle_tory_z_odstepem_nm(nm):
    a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    offset_deg = (nm * 1852.0) / 111320.0
    b = _prosty_tor(50.0 + offset_deg, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    return a, b


def test_domyslny_airspace_to_en_route_5nm_1000ft(timdr):
    a, b = _rownolegle_tory_z_odstepem_nm(4.0)
    domyslny = timdr.conflict_alert(a, b, lookahead_seconds=60.0)
    jawny_en_route = timdr.conflict_alert(a, b, airspace="en_route", lookahead_seconds=60.0)
    assert domyslny["horizontal_threshold_nm"] == jawny_en_route["horizontal_threshold_nm"] == 5.0
    assert domyslny["vertical_threshold_ft"] == jawny_en_route["vertical_threshold_ft"] == 1000.0
    assert domyslny["airspace"] == "en_route"


def test_preset_tma_ma_mniejszy_prog_niz_en_route(timdr):
    """Odstep 4 NM: konflikt wg en-route (prog 5NM), bezpiecznie wg TMA
    (prog 3NM) - to samo geometrycznie, rozne wnioski w zaleznosci od
    przestrzeni, dokladnie jak w realnym ATC."""
    a, b = _rownolegle_tory_z_odstepem_nm(4.0)
    en_route = timdr.conflict_alert(a, b, airspace="en_route", lookahead_seconds=60.0)
    tma = timdr.conflict_alert(a, b, airspace="tma", lookahead_seconds=60.0)
    assert en_route["conflict"] is True
    assert tma["conflict"] is False
    assert tma["horizontal_threshold_nm"] == 3.0


def test_jawne_progi_nadpisuja_preset(timdr):
    a, b = _rownolegle_tory_z_odstepem_nm(4.0)
    wynik = timdr.conflict_alert(a, b, airspace="tma", horizontal_nm=10.0, lookahead_seconds=60.0)
    assert wynik["horizontal_threshold_nm"] == 10.0  # jawna wartosc, nie 3.0 z presetu
    assert wynik["conflict"] is True  # 4NM < 10NM


def test_nieznany_airspace_rzuca_wyjatek(timdr):
    a, b = _rownolegle_tory_z_odstepem_nm(4.0)
    with pytest.raises(ValueError):
        timdr.conflict_alert(a, b, airspace="strefa_ktora_nie_istnieje")


def test_fleet_wykrywa_wlasciwe_pary(timdr):
    """3 tory: A i B na kursie kolizyjnym (konflikt), C rownolegle daleko
    od obu (bezpieczny) - fleet powinien zglosic dokladnie jedna pare."""
    a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200)
    c = _prosty_tor(50.0 + 18520.0 / 111320.0 * 5, 19.5, 10000.0, heading_deg=90, speed_mps=200)  # 50 NM dalej

    wyniki = timdr.conflict_alert_fleet({"A": a, "B": b, "C": c}, lookahead_seconds=180.0)

    assert len(wyniki) == 1
    para = {wyniki[0]["a_id"], wyniki[0]["b_id"]}
    assert para == {"A", "B"}
    assert wyniki[0]["conflict"] is True


def test_fleet_include_safe_zwraca_wszystkie_pary(timdr):
    a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200)
    c = _prosty_tor(50.0 + 18520.0 / 111320.0 * 5, 19.5, 10000.0, heading_deg=90, speed_mps=200)

    wyniki = timdr.conflict_alert_fleet({"A": a, "B": b, "C": c}, lookahead_seconds=180.0, include_safe=True)
    assert len(wyniki) == 3  # C(3,2) = 3 pary
    liczba_konfliktow = sum(1 for w in wyniki if w["conflict"])
    assert liczba_konfliktow == 1


def test_fleet_posortowany_po_pilnosci(timdr):
    """Para z krotszym time_to_conflict_s powinna byc pierwsza."""
    a = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200)  # bliski konflikt
    # druga para: taki sam kurs kolizyjny, ale dalej od siebie (pozniejszy konflikt)
    c = _prosty_tor(51.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    d = _prosty_tor(51.0, 22.5, 10000.0, heading_deg=270, speed_mps=200)

    wyniki = timdr.conflict_alert_fleet(
        {"A": a, "B": b, "C": c, "D": d}, lookahead_seconds=600.0
    )
    czasy = [w["time_to_conflict_s"] for w in wyniki if w["conflict"]]
    assert czasy == sorted(czasy)
