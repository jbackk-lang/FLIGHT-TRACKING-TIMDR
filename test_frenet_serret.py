"""
test_frenet_serret.py -- testy dla frenet_serret()/twist_3d() (kappa/tau
prawdziwej trajektorii 3D), dodanych po sesji projektowej w której
zweryfikowano dokladnie te same wzory (kappa, tau z r', r'', r''') na
helisie Boerdijka-Coxetera i na czystej analitycznej helisie kolowej
(zob. skill timdr-signal-framework SS16). Tutaj to samo, ale jako formalny
test regresyjny na module produkcyjnym, przechodzacym przez pelny pipeline
geo-projekcji (lat/lon -> lokalna plaszczyzna styczna), nie tylko na
surowych wspolrzednych ENU.

Kluczowe, uczciwie udokumentowane ograniczenie: tau wymaga trzeciej
pochodnej pozycji (jerk). Surowe roznicowanie (np.gradient x3) okazalo sie
w praktyce bezuzyteczne pod realistycznym szumem GPS (std tau ~40x wieksze
od sygnalu przy szumie 0.5m) - stad wygladzanie Savitzky-Golay z oknem
dobieranym w SEKUNDACH wewnatrz frenet_serret(). Test
test_tau_degraduje_pod_szumem dokumentuje ten kompromis wprost, zamiast
udawac ze go nie ma.
"""
import numpy as np
import pytest
from timdr_flight import TIMDRFlight


@pytest.fixture
def timdr():
    return TIMDRFlight()


def _helisa_track(timdr, R=3000.0, omega=0.05, c=15.0, n=400, T=120.0,
                   lat0=50.0, lon0=20.0, alt0=1000.0, noise_m=0.0, seed=0):
    t = np.linspace(0, T, n)
    x = R * np.cos(omega * t)
    y = R * np.sin(omega * t)
    z = alt0 + c * t
    if noise_m > 0:
        rng = np.random.default_rng(seed)
        x = x + rng.normal(0, noise_m, n)
        y = y + rng.normal(0, noise_m, n)
        z = z + rng.normal(0, noise_m, n)
    lat, lon = timdr._inverse_local_xy(np.column_stack([x, y]), lat0, lat0, lon0)
    return np.column_stack([lat, lon, z, t]), (R, omega, c)


def test_frenet_serret_helisa_kappa_tau_dokladne(timdr):
    """Waliduje kappa/tau na trajektorii o znanym analitycznie rozwiazaniu
    (helisa kolowa), przechodzacej przez PELNY pipeline (lat/lon -> ENU),
    nie tylko surowe wspolrzedne."""
    track, (R, omega, c) = _helisa_track(timdr)
    kappa, tau = timdr.frenet_serret(track, window_seconds=6.0, poly=3)

    kappa_true = R * omega**2 / (R**2 * omega**2 + c**2)
    tau_true = c * omega / (R**2 * omega**2 + c**2)

    mid = slice(len(track) // 4, -len(track) // 4)
    kappa_err = abs(np.median(kappa[mid]) - kappa_true) / kappa_true
    tau_err = abs(np.median(tau[mid]) - tau_true) / tau_true

    assert kappa_err < 0.05, f"blad kappa {kappa_err:.3%} (kappa={np.median(kappa[mid])}, prawda={kappa_true})"
    assert tau_err < 0.05, f"blad tau {tau_err:.3%} (tau={np.median(tau[mid])}, prawda={tau_true})"


def test_tau_degraduje_pod_szumem_ale_nie_eksploduje(timdr):
    """Uczciwy test ograniczenia: przy realistycznym szumie GPS (0.5m)
    blad tau rosnie, ale ze zgladzaniem Savitzky-Golay (okno w sekundach)
    pozostaje rzedu wielkosci sygnalu, a nie 10-40x wiekszy jak przy
    surowym trzykrotnym roznicowaniu (zmierzone w sesji projektowej)."""
    track, (R, omega, c) = _helisa_track(timdr, noise_m=0.5, seed=1)
    kappa, tau = timdr.frenet_serret(track, window_seconds=6.0, poly=3)
    tau_true = c * omega / (R**2 * omega**2 + c**2)

    mid = slice(len(track) // 4, -len(track) // 4)
    tau_std = np.std(tau[mid])
    # std powinno pozostac tego samego rzedu wielkosci co sam sygnal,
    # nie 10-40x wieksze (to bylby dowod ze wygladzanie nie dziala)
    assert tau_std < 10 * abs(tau_true), (
        f"tau_std={tau_std} zbyt duze wzgledem sygnalu tau={tau_true} - "
        f"wygladzanie nie tlumi szumu jak oczekiwano"
    )


def test_twist_3d_lot_prostoliniowy_bez_falszywych_alarmow(timdr):
    """Lot w linii prostej i na stalej wysokosci: tau powinno byc bliskie
    zeru wszedzie, wiec adaptacyjny prog (wzgledny) nie powinien flagowac
    nic - podloga (floor_frac) ma temu zapobiegac."""
    n = 120
    t = np.linspace(0, 600, n)
    lat = np.full(n, 50.0) + 0.0005 * t / 600.0  # powolny, plynny ruch na polnoc
    lon = np.full(n, 20.0) + 0.0005 * t / 600.0
    alt = np.full(n, 10000.0)
    track = np.column_stack([lat, lon, alt, t])

    idx = timdr.twist_3d(track, tau_factor=3.0)
    assert len(idx) == 0, f"falszywe alarmy skretu na prostym locie: {idx}"


def test_twist_3d_wykrywa_ostry_manewr(timdr):
    """Do plynnego, lagodnego zakretu (staly promien, stala torsja)
    wstrzykujemy krotki, gwaltowny 'korkociag' (duza, szybka zmiana
    kierunku WEKTORA normalnego, nie tylko kursu) i sprawdzamy, czy
    twist_3d go lapie, a caly reszta lotu (plynny zakret) zostaje
    niezaflagowana."""
    R, omega, c = 4000.0, 0.04, 5.0
    n_pre, n_maneuver, n_post = 80, 12, 80
    t = np.linspace(0, (n_pre + n_maneuver + n_post - 1) * 1.5, n_pre + n_maneuver + n_post)

    x = np.zeros(len(t)); y = np.zeros(len(t)); z = np.zeros(len(t))
    x[:n_pre] = R * np.cos(omega * t[:n_pre])
    y[:n_pre] = R * np.sin(omega * t[:n_pre])
    z[:n_pre] = 8000 + c * t[:n_pre]

    # gwaltowny manewr: szybka oscylacja w osi pionowej (barrel-roll-like)
    # nalozona na kontynuacje lagodnego zakretu, zeby nie wprowadzac
    # sztucznego skoku pozycji (ta sama zasada co przy testach recovery
    # w innych repo TIMDR - kontynuuj proces bazowy, nie zaczynaj od nowa)
    seg = slice(n_pre, n_pre + n_maneuver)
    t_seg = t[seg]
    x[seg] = R * np.cos(omega * t_seg)
    y[seg] = R * np.sin(omega * t_seg)
    z[seg] = 8000 + c * t_seg + 400 * np.sin(2 * np.pi * (t_seg - t[n_pre]) / (1.5 * (n_maneuver - 1)))

    t_post0 = t[n_pre + n_maneuver]
    x[n_pre + n_maneuver:] = R * np.cos(omega * t[n_pre + n_maneuver:])
    y[n_pre + n_maneuver:] = R * np.sin(omega * t[n_pre + n_maneuver:])
    z[n_pre + n_maneuver:] = z[n_pre + n_maneuver - 1] + c * (t[n_pre + n_maneuver:] - t[n_pre + n_maneuver - 1])

    lat, lon = timdr._inverse_local_xy(np.column_stack([x, y]), 50.0, 50.0, 20.0)
    track = np.column_stack([lat, lon, z, t])

    idx = timdr.twist_3d(track, tau_factor=3.0, window_seconds=6.0)
    assert len(idx) > 0, "manewr nie zostal wykryty"

    # wiekszosc flag powinna wypasc w oknie manewru (+- margines na
    # wygladzanie), nie rozrzucona losowo po calym plynnym zakrecie
    in_maneuver = np.sum((idx >= n_pre - 5) & (idx <= n_pre + n_maneuver + 5))
    assert in_maneuver / len(idx) > 0.5, (
        f"flagi rozrzucone poza oknem manewru: {idx} (manewr: {n_pre}-{n_pre+n_maneuver})"
    )
