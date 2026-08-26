import numpy as np
import pytest
from timdr_flight import TIMDRFlight


@pytest.fixture
def timdr():
    return TIMDRFlight()


def test_walidacja_ksztaltu(timdr):
    with pytest.raises(ValueError):
        timdr.twist([[50.0, 19.9, 1000], [50.1, 20.0, 1100]])  # brak kolumny t


def test_walidacja_min_punktow(timdr):
    with pytest.raises(ValueError):
        timdr.twist([[50.0, 19.9, 1000, 0]])


def test_walidacja_nierosnacy_czas(timdr):
    track = [[50.0, 19.9, 1000, 0], [50.1, 20.0, 1100, 0]]  # dt = 0
    with pytest.raises(ValueError):
        timdr.twist(track)


def test_walidacja_antimeridian(timdr):
    track = [[10.0, 179.9, 1000, 0], [10.0, -179.9, 1000, 10]]  # skok 359.8 st
    with pytest.raises(ValueError):
        timdr.twist(track)


def test_kurs_naiwny_kontra_poprawny_roznia_sie(timdr):
    """
    Regresja dla bugu jednostek geograficznych: surowe stopnie lat/lon nie
    sa izotropowe. Sprawdzamy, ze zaimplementowany kurs (z korekta
    cos(lat)) NIE jest rownowazny naiwnemu arctan2 na surowych stopniach
    dla lotu na srednich szerokosciach geograficznych.
    """
    track = [
        [50.0, 19.9, 1000, 0],
        [50.01, 19.91, 1200, 10],
        [50.03, 19.93, 1500, 20],
        [50.06, 19.96, 1800, 30],
        [50.10, 20.00, 2000, 40],
    ]
    diag = timdr.diagnostics(track)
    course = diag["course_deg"]

    tr = np.array(track)
    v_naive = np.gradient(tr[:, :2], axis=0)
    course_naive = np.rad2deg(np.arctan2(v_naive[:, 1], v_naive[:, 0])) % 360.0

    # kursy powinny sie realnie roznic (o >5 st) - to dowod ze poprawka
    # cos(lat) faktycznie cos zmienia, a nie jest no-opem
    assert np.any(np.abs((course - course_naive + 180) % 360 - 180) > 5)


def test_lot_na_poludnie_bez_falszywego_twistu(timdr):
    """
    Regresja dla bugu zawijania kata przy kursie ~180 st (lot na
    poludnie). Rzeczywista zmiana kursu ~0.1-0.6 st/krok, prog 20 st -
    nie powinno byc zadnych alarmow.
    """
    track = [
        [50.00, 20.000, 1000, 0],
        [49.90, 20.0005, 1000, 10],
        [49.80, 19.9995, 1000, 20],
        [49.70, 20.0005, 1000, 30],
        [49.60, 19.9995, 1000, 40],
    ]
    result = timdr.twist(track)
    assert len(result["direction_twist"]) == 0, (
        f"falszywe alarmy kursu z powodu zawijania kata: {result['direction_twist']}"
    )


def test_ostry_zwrot_wykryty(timdr):
    # lot na wschod, potem ostry zwrot 90 st na poludnie
    track = [
        [50.0, 19.0, 1000, 0],
        [50.0, 19.2, 1000, 10],
        [50.0, 19.4, 1000, 20],
        [49.8, 19.4, 1000, 30],
        [49.6, 19.4, 1000, 40],
        [49.4, 19.4, 1000, 50],
    ]
    result = timdr.twist(track)
    assert len(result["direction_twist"]) > 0


def test_prog_wysokosci_niezalezny_od_probkowania(timdr):
    """
    Regresja dla bugu jednostek predkosci pionowej: ta sama fizyczna
    predkosc wznoszenia (~20 m/s) wykryta niezaleznie od interwalu
    probkowania (10s vs 40s).
    """
    fast_sample = [
        [50.0, 20.0, 1000, 0],
        [50.0, 20.0, 1200, 10],
        [50.0, 20.0, 1400, 20],
    ]
    slow_sample = [
        [50.0, 20.0, 1000, 0],
        [50.0, 20.0, 1800, 40],
        [50.0, 20.0, 2600, 80],
    ]
    r1 = timdr.twist(fast_sample, climb_rate_thresh_mps=15.0)
    r2 = timdr.twist(slow_sample, climb_rate_thresh_mps=15.0)
    # obie trasy maja predkosc wznoszenia 20 m/s -> obie powinny przekroczyc prog
    assert len(r1["altitude_twist"]) > 0
    assert len(r2["altitude_twist"]) > 0


def test_predict_ruch_jednostajny(timdr):
    # lot ze stala predkoscia pozioma, bez wznoszenia -> predykcja powinna
    # kontynuowac trase w przyblizeniu liniowo
    track = [
        [50.0, 19.0 + 0.01 * i, 1000, float(i * 10)] for i in range(5)
    ]
    pred = timdr.predict(track, steps=3)
    assert pred.shape == (3, 3)
    assert not np.any(np.isnan(pred))
    # lon powinna nadal rosnac w przyblizeniu tym samym krokiem
    assert pred[0, 1] > track[-1][1]
    assert pred[1, 1] > pred[0, 1]


def test_trm_reduce_zachowuje_konce(timdr):
    track = [
        [50.0, 20.0, 1000, 0],
        [50.1, 20.1, 1100, 10],
        [50.0, 20.0, 1000, 20],
        [50.1, 20.1, 1100, 30],
        [50.0, 20.0, 1000, 40],
    ]
    smooth = timdr.trm_reduce(track)
    assert np.allclose(smooth[0], track[0][:3])
    assert np.allclose(smooth[-1], track[-1][:3])
    assert not np.any(np.isnan(smooth))


def test_diagnostics_jednostki_sensowne(timdr):
    track = [
        [50.0, 19.9, 1000, 0],
        [50.01, 19.91, 1200, 10],
        [50.03, 19.93, 1500, 20],
        [50.06, 19.96, 1800, 30],
        [50.10, 20.00, 2000, 40],
    ]
    diag = timdr.diagnostics(track)
    assert np.all(diag["course_deg"] >= 0) and np.all(diag["course_deg"] < 360)
    assert np.all(diag["ground_speed_mps"] >= 0)
    # predkosc wznoszenia rzedu 20-30 m/s -> kilkaset-kilka tysiecy ft/min
    assert np.all(np.abs(diag["climb_rate_fpm"]) < 10000)
