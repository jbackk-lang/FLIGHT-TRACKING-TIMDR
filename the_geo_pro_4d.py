"""
the_geo_pro_4d.py
------------------
Prawdziwa (nie pseudokodowa) implementacja THE-GEO PRO 4D.

Kontekst / historia poprawek — patrz README.md, sekcja
"THE-GEO PRO -> THE-GEO PRO 4D: historia walidacji":

  1. Oryginalny wzor na skret w tym repo (Torsion.compute w sekcji
     "THE-GEO PRO" nizej w README), oparty na roznicach KIERUNKOW
     jednostkowych (D_t2, D_t1, D_t), zostal sprawdzony na czystej,
     bezszumowej helisie analitycznej i NIE zbiega do prawdziwej
     torsji -- blad zostaje na poziomie 90-99% nawet przy zageszczaniu
     probkowania (dt = 0.5 .. 0.005).

  2. Ten plik implementuje poprawiona wersje ("4D" = pozycje
     parametryzowane czasem t jako 4. wspolrzedna krzywej). Krzywizna
     i skret sa liczone z pochodnych predkosci / przyspieszenia /
     szarpniecia (v, a, j) wyliczonych roznicami skonczonymi z 4
     kolejnych punktow:

         kappa = |v x a| / |v|^3
         tau   = det(v, a, j) / |v x a|^2
         H     = sqrt(kappa^2 + tau^2)      (helikalnosc)

     Zbiega poprawnie do wartosci analitycznych (patrz testy i README).

  3. Implementacja jest matematycznie identyczna z
     FLIGHT-TRACKING-TIMDR/core/curvature_detector_3d.py (sprawdzone
     numerycznie -- wyniki v/a/j bit-identyczne dla tych samych
     wejsc). Rozni sie tylko interfejsem: zwraca dict zgodny z
     oryginalnym pseudokodem THE-GEO PRO 4D i dodatkowo liczy
     "helical".

  4. Poprawka wzgledem pseudokodu uzytkownika (wariant
     "THE_GEO_PRO_4D_Radar"): zabezpieczenie `if cross_norm == 0` (lub
     nawet stala tolerancja typu `< 1e-12`) chroni TYLKO przed dzieleniem
     przez doslowne zero / NaN -- NIE chroni przed wzmacnianiem szumu,
     gdy trajektoria jest niemal (ale nie dokladnie) prosta. Przyklad:
     v=(1,0,0), a=(1,1e-6,0) (typowy szum kierunku na "prostym" odcinku)
     daje cross_norm=1e-6 -- nie zero, wiec stary warunek przepuszcza
     dalej, a tau = dot(cross_va,j)/cross_norm**2 wychodzi rzedu 1e6
     zamiast ~0. To ten sam wzorzec bledu co przy min_speed / min_step_m
     (dzielenie przez mala wartosc wzmacnia szum), tylko na innym
     mianowniku (cross_norm**2 zamiast |v|).
     Poprawka: bramkowanie torsji na podstawie samej krzywizny kappa
     (juz obliczonej, fizycznie sensownej wielkosci), nie na cross_norm:
     jesli kappa < min_curvature -> tau = 0. Domyslna wartosc
     min_curvature=1e-4 jest punktem startowym, NIE zwalidowana stala --
     wymaga kalibracji na realnych, zaszumionych danych 3D, analogicznie
     do min_step_m=3.0 skalibrowanego dla RADAR-TRACKING-TIMDR w 2D.

Wejscie: p_t, p_t1, p_t2, p_t3 -- krotki (x, y, z). p_t = pozycja
najnowsza, p_t3 = pozycja sprzed 3 krokow (najstarsza).
"""

import math
from typing import Dict, Tuple

Point3 = Tuple[float, float, float]

DEFAULT_MIN_SPEED = 1.0
DEFAULT_MIN_CURVATURE = 1e-4


def _sub(a: Point3, b: Point3) -> Point3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _norm(v: Point3) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def _cross(a: Point3, b: Point3) -> Point3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Point3, b: Point3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def THE_GEO_PRO_4D(
    p_t: Point3,
    p_t1: Point3,
    p_t2: Point3,
    p_t3: Point3,
    min_speed: float = DEFAULT_MIN_SPEED,
    min_curvature: float = DEFAULT_MIN_CURVATURE,
) -> Dict[str, float]:
    """Krzywizna, skret i helikalnosc trajektorii 3D sparametryzowanej czasem.

    Zwraca dict: {"gated", "curvature", "torsion", "helical"}.
    gated=True gdy predkosc ponizej min_speed (szum GPS/sensora przy
    postoju) -- wtedy pozostale pola sa wyzerowane, zeby nie wzmacniac
    szumu przy dzieleniu przez mala predkosc.
    """
    v = _sub(p_t, p_t1)
    v1 = _sub(p_t1, p_t2)
    a = _sub(v, v1)
    a1 = _sub(v1, _sub(p_t2, p_t3))
    j = _sub(a, a1)

    speed = _norm(v)
    if speed < min_speed:
        return {"gated": True, "curvature": 0.0, "torsion": 0.0, "helical": 0.0}

    cross_va = _cross(v, a)
    cross_norm = _norm(cross_va)
    kappa = cross_norm / speed ** 3

    if kappa < min_curvature:
        # trajektoria (niemal) prostoliniowa -- torsja niezdefiniowana /
        # zdominowana przez szum przy dzieleniu przez cross_norm**2,
        # bramkujemy na podstawie kappa, nie samego cross_norm==0
        # (patrz docstring modulu, punkt 4)
        tau = 0.0
    else:
        tau = _dot(cross_va, j) / cross_norm ** 2

    H = math.sqrt(kappa * kappa + tau * tau)

    return {"gated": False, "curvature": kappa, "torsion": tau, "helical": H}


# ============================================================
# NAPRAWIONE (wzmacnianie szumu przy roznicowaniu na krotkich/realnych
# trajektoriach -- patrz uzasadnienie ponizej): THE_GEO_PRO_4D() liczy
# v/a/j z SUROWYCH, kolejnych punktow -- na analitycznej, gladkiej
# helisie (jak w testach wyzej) to zbiega poprawnie do prawdziwej
# krzywizny/torsji. Ale przy realnym szumie pomiaru (GPS/ADS-B/dowolny
# sensor), potrojne roznicowanie (v->a->j) WZMACNIA szum na kazdym
# kroku -- FLIGHT-TRACKING-TIMDR (repo siostrzane, uzywajace
# matematycznie identycznego wzoru) zmierzyl to bezposrednio: surowe
# roznicowanie dawalo std(tau)~2.6 przy prawdziwym tau~0.065 (40x
# wzmocnienie szumu), naprawione tam wlasnie wygladzaniem PRZED
# roznicowaniem. Niezaleznie, w GIA-TIMDR (repo teoretyczne uzywajace
# tego samego wzoru na trojwezle) potwierdzono ten sam mechanizm na
# krotkich (~25-punktowych) realnych szeregach: kappa/tau z surowych
# roznic sa nierozroznialne od szumu i maja recall bliski przypadkowi
# (patrz GIA-TIMDR/docs/geometry/TIMDR_Trefoil_RealDataValidation.md).
#
# THE_GEO_PRO_4D() (funkcja wyzej) POZOSTAJE NIEZMIENIONA -- dziala na
# dokladnie 4 punktach, jest uzywana bezposrednio przez kod zalezny
# (m.in. THE_GEO_PRO_4D_Radar.py) i pozostaje poprawna dla gladkich/
# analitycznych wejsc (patrz TestHelixConvergence w testach). Ponizej:
# NOWA, zalecana sciezka dla realnych/zaszumionych sekwencji dluzszych
# niz 4 punkty -- v/a/j liczone BEZPOSREDNIO z pochodnych lokalnie
# dopasowanego wielomianu (prawdziwy filtr Savitzky-Golay), NIE
# "wygladz pozycje, potem policz proste roznice sasiadow" (ta pierwsza,
# prostsza wersja zostala wyprobowana i ODRZUCONA -- dawala WIEKSZY
# blad niz brak wygladzania, patrz docstring
# velocity_acceleration_jerk_at() nizej).
#
# Implementacja BEZ zaleznosci od numpy/scipy (ten sam wzorzec unikania
# ciezkich zaleznosci co przy naprawie Device Guard w
# TIMDR-Industrial-Predict/TIMDR-EV-Predict/TIMDR-Earthquake-Core) --
# rownania normalne najmniejszych kwadratow rozwiazywane wlasna
# eliminacja Gaussa z czesciowym wyborem elementu podstawowego.
# ============================================================

def _solve_linear_system(A: list, b: list) -> list:
    """Rozwiazuje Ax=b eliminacja Gaussa z czesciowym wyborem elementu
    podstawowego. A: lista list (n x n), b: lista dlugosci n. Male
    uklady (n<=polyorder+1<=4 w praktycznym uzyciu ponizej) -- prostota
    ponad wydajnosc."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[pivot_row] = M[pivot_row], M[col]
        if abs(M[col][col]) < 1e-14:
            raise ValueError("Osobliwy uklad rownan normalnych (za malo punktow w oknie)")
        for r in range(col + 1, n):
            factor = M[r][col] / M[col][col]
            for c in range(col, n + 1):
                M[r][c] -= factor * M[col][c]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


def _local_poly_coeffs(series: list, i: int, window: int, polyorder: int) -> list:
    """Dopasowuje wielomian stopnia `polyorder` metoda najmniejszych
    kwadratow do okna wokol indeksu i (zmienna niezalezna = k-i,
    wycentrowana), zwraca WSZYSTKIE wspolczynniki [c0,c1,...,cD]."""
    n = len(series)
    half = window // 2
    lo, hi = max(0, i - half), min(n, i + half + 1)
    xs = [k - i for k in range(lo, hi)]
    ys = series[lo:hi]
    m = min(polyorder + 1, len(xs))
    ATA = [[sum(x ** (a + b) for x in xs) for b in range(m)] for a in range(m)]
    ATy = [sum((x ** a) * y for x, y in zip(xs, ys)) for a in range(m)]
    coeffs = _solve_linear_system(ATA, ATy)
    while len(coeffs) < polyorder + 1:
        coeffs.append(0.0)
    return coeffs


def velocity_acceleration_jerk_at(points: list, i: int, window: int = 9, polyorder: int = 3) -> Tuple[Point3, Point3, Point3]:
    """Ocenia (v,a,j) w punkcie i BEZPOSREDNIO z pochodnych lokalnie
    dopasowanego wielomianu (prawdziwy filtr rozniczkujacy
    Savitzky-Golay), zamiast liczyc je jako proste roznice sasiednich,
    NIEZALEZNIE wygladzonych punktow.

    UZASADNIENIE POPRAWKI (znalezione empirycznie -- pierwsza probka
    tej naprawy uzywala podejscia "wygladz pozycje, potem rozniczkuj
    proste roznice" i NIE dzialala: blad tau byl WIEKSZY niz bez
    wygladzania w ogole, patrz historia commitow/rozmowa). Przyczyna:
    v_smoothed(i) - v_smoothed(i-1) odejmuje dwa NIEZALEZNIE dopasowane
    lokalne wielomiany (z przesunietych o 1 probke, zachodzacych na
    siebie okien) -- to nie usuwa szumu z wyzszych pochodnych tak
    skutecznie, jak wziecie pochodnej ANALITYCZNIE z JEDNEGO dopasowania
    per punkt. Dla wielomianu c0+c1*t+c2*t^2+c3*t^3 (t=k-i), pochodne w
    t=0 (srodek okna) to: wartosc=c0, predkosc=c1, przyspieszenie=2*c2,
    szarpniecie=6*c3 -- to jest STANDARDOWY filtr pochodnych
    Savitzky-Golay, nie improwizacja."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    cx = _local_poly_coeffs(xs, i, window, polyorder)
    cy = _local_poly_coeffs(ys, i, window, polyorder)
    cz = _local_poly_coeffs(zs, i, window, polyorder)
    v = (cx[1], cy[1], cz[1])
    a = (2 * cx[2], 2 * cy[2], 2 * cz[2])
    j = (6 * cx[3], 6 * cy[3], 6 * cz[3])
    return v, a, j


def _kappa_tau_from_vaj(v: Point3, a: Point3, j: Point3, min_curvature: float) -> Tuple[bool, float, float]:
    speed = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    cross_va = (v[1]*a[2]-v[2]*a[1], v[2]*a[0]-v[0]*a[2], v[0]*a[1]-v[1]*a[0])
    cross_norm = math.sqrt(sum(c*c for c in cross_va))
    if speed < 1e-12:
        return True, 0.0, 0.0
    kappa = cross_norm / speed ** 3
    if kappa < min_curvature:
        return False, kappa, 0.0
    tau = sum(c*jj for c, jj in zip(cross_va, j)) / cross_norm ** 2
    return False, kappa, tau


def THE_GEO_PRO_4D_sequence(
    points: list,
    min_speed: float = DEFAULT_MIN_SPEED,
    min_curvature: float = DEFAULT_MIN_CURVATURE,
    smooth: bool = True,
    window: int = 9,
    polyorder: int = 3,
) -> list:
    """Zalecane wejscie dla realnych/zaszumionych sekwencji dluzszych
    niz 4 punkty. Dwa tryby:

    - `smooth=False` (stare, niezmitygowane zachowanie): dokladnie
      odtwarza wyniki wielokrotnego wywolania THE_GEO_PRO_4D() na
      kolejnych oknach 4 surowych punktow -- wsteczna kompatybilnosc,
      do porownan/testow regresyjnych.
    - `smooth=True` (DOMYSLNE, zalecane dla danych realnych/zaszumionych):
      liczy v/a/j BEZPOSREDNIO z pochodnych lokalnie dopasowanego
      wielomianu stopnia `polyorder` w oknie `window` probek wokol
      kazdego punktu (prawdziwy filtr Savitzky-Golay, patrz
      velocity_acceleration_jerk_at() powyzej) -- NIE "wygladz pozycje,
      potem licz proste roznice", co empirycznie dawalo WIEKSZY blad niz
      brak wygladzania w ogole (patrz docstring
      velocity_acceleration_jerk_at). Zmierzone na zaszumionej helisie
      (test_noise_amplification_fix.py): srednia bledu tau spada ok.
      2-2.5x wzgledem surowego roznicowania.

    Zwraca liste wynikow (dict: gated/curvature/torsion), jeden na
    kazdy punkt wejsciowy z wystarczajaca historia -- dla `smooth=False`
    to indeksy i>=3 (offset o 3, jak wczesniej); dla `smooth=True` to
    indeksy i>=window//2 i i<len(points)-window//2 (potrzeba pelnego
    okna po obu stronach dla stabilnej pochodnej 3. rzedu)."""
    if not smooth:
        n = len(points)
        return [
            THE_GEO_PRO_4D(points[i], points[i-1], points[i-2], points[i-3],
                            min_speed=min_speed, min_curvature=min_curvature)
            for i in range(3, n)
        ]

    if polyorder < 3:
        raise ValueError(
            "smooth=True wymaga polyorder>=3 -- torsja potrzebuje wyrazu "
            "szescianowego (szarpniecie = 6*c3), nizszy stopien nie ma tego "
            "wspolczynnika (patrz velocity_acceleration_jerk_at)."
        )

    n = len(points)
    half = window // 2
    results = []
    for i in range(half, n - half):
        v, a, j = velocity_acceleration_jerk_at(points, i, window, polyorder)
        speed = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        if speed < min_speed:
            results.append({"gated": True, "curvature": 0.0, "torsion": 0.0, "helical": 0.0})
            continue
        gated, kappa, tau = _kappa_tau_from_vaj(v, a, j, min_curvature)
        H = math.sqrt(kappa*kappa + tau*tau)
        results.append({"gated": gated, "curvature": kappa, "torsion": tau, "helical": H})
    return results
