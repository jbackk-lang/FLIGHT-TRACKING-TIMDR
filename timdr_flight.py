"""
FLIGHT-TRACKING-TIMDR / timdr_flight.py
========================================
Moduł analizy toru lotu: gradient ruchu (TIMDR-flow), detekcja nagłych
zmian kursu i wysokości (twist), redukcja szumu toru (TRM) oraz prosta
predykcja trajektorii.

Wejście: tor lotu jako lista/tablica punktów [lat, lon, alt, t]
  - lat, lon: stopnie dziesiętne (WGS84)
  - alt: metry (jeśli masz dane w stopach, przelicz przed wywołaniem: m = ft * 0.3048)
  - t: sekundy (musi ściśle rosnąć)
"""

import numpy as np

EARTH_RADIUS_M = 6371000.0  # średni promień Ziemi - wystarczający dla torów
                             # regionalnych; patrz ograniczenia w __doc__ klasy


class TIMDRFlight:
    """
    Ograniczenia geodezyjne (ważne, przeczytaj przed użyciem na realnych
    danych ADS-B):

    - Pozycja jest rzutowana na lokalną płaszczyznę styczną (przybliżenie
      equirectangular) względem średniej szerokości geograficznej toru.
      To przybliżenie jest dobre dla torów lokalnych/regionalnych
      (rzędu do kilkuset km). Dla lotów transoceanicznych/długodystansowych
      błąd rzutowania rośnie i lepiej dzielić trasę na segmenty albo użyć
      właściwej projekcji geodezyjnej (np. pyproj).
    - Tor przecinający południk 180° (antimeridian) NIE jest obsługiwany
      poprawnie - `_validate` zgłosi wyjątek, jeśli wykryje skok długości
      geograficznej > 180° między kolejnymi punktami, zamiast cicho zwracać
      błędny wynik.
    """

    def __init__(self):
        pass

    # ------------------------------------------------------------
    # walidacja i geometria wspólna
    # ------------------------------------------------------------
    @staticmethod
    def _validate(track, min_points=2):
        tr = np.asarray(track, dtype=np.float64)
        if tr.ndim != 2 or tr.shape[1] != 4:
            raise ValueError(
                f"track musi mieć kształt (N, 4) [lat, lon, alt, t], dostano {tr.shape}"
            )
        if len(tr) < min_points:
            raise ValueError(
                f"track musi mieć co najmniej {min_points} punkty, dostano {len(tr)}"
            )
        t = tr[:, 3]
        if np.any(np.diff(t) <= 0):
            # POPRAWKA: tak jak w TIMDR-Radar-Module - dt <= 0 daje dzielenie
            # przez zero w np.gradient(pos, t, axis=0).
            raise ValueError(
                "znaczniki czasu (kolumna t) muszą być ściśle rosnące (dt > 0)"
            )
        lon = tr[:, 1]
        if np.any(np.abs(np.diff(lon)) > 180.0):
            # Tor przecinający +-180 dlugosci geograficznej (antimeridian)
            # daje bez korekty ogromny, falszywy skok pozycji/predkosci.
            raise ValueError(
                "wykryto skok długości geograficznej > 180° między kolejnymi "
                "punktami - prawdopodobnie tor przecina antimeridian (±180°); "
                "ten moduł tego nie obsługuje, trzeba podzielić tor na segmenty"
            )
        return tr

    @staticmethod
    def _project_local_xy(lat, lon, lat0):
        """
        Equirectangular: lokalna plaszczyzna styczna w metrach.
        x = wschod (dodatni), y = polnoc (dodatni).

        POPRAWKA (bug jednostek geograficznych): oryginalny kod liczył
        gradient bezpośrednio na surowych stopniach [lat, lon], traktując
        1 stopień lat i 1 stopień lon jako tę samą "odległość". To
        nieprawda poza równikiem: 1° długości geograficznej odpowiada
        ok. 111.32 km * cos(lat) na powierzchni Ziemi, czyli maleje wraz
        ze wzrostem szerokości. Zweryfikowano empirycznie na torze z
        przykładu (lat ok. 50°N): kurs liczony z surowych stopni dawał
        45.0°, podczas gdy prawdziwy kurs (po korekcie cos(lat)) to
        57.29° - błąd 12.3°, zdecydowanie za duży dla zastosowań ATC.
        """
        x = np.deg2rad(lon - lon[0]) * EARTH_RADIUS_M * np.cos(np.deg2rad(lat0))
        y = np.deg2rad(lat - lat[0]) * EARTH_RADIUS_M
        return np.column_stack([x, y])

    @staticmethod
    def _inverse_local_xy(xy, lat0, lat_ref, lon_ref):
        lat = lat_ref + np.rad2deg(xy[:, 1] / EARTH_RADIUS_M)
        lon = lon_ref + np.rad2deg(xy[:, 0] / (EARTH_RADIUS_M * np.cos(np.deg2rad(lat0))))
        return lat, lon

    # --- 1. TIMDR-flow: gradient ruchu lotu ---
    def timdr_flow(self, track):
        """
        track: lista punktów [lat, lon, alt, t]
        zwraca: TIMDR-flow (gradient prędkości+przyspieszenia) w metrach,
        jako [flow_x_wschod, flow_y_polnoc, flow_alt] - spójne jednostki
        (m/s^2 na krok czasu) dla wszystkich trzech składowych.
        """
        tr = self._validate(track)
        lat, lon, alt, t = tr[:, 0], tr[:, 1], tr[:, 2], tr[:, 3]
        lat0 = float(np.mean(lat))

        xy = self._project_local_xy(lat, lon, lat0)
        pos_m = np.column_stack([xy, alt])  # metry we wszystkich 3 osiach

        v = np.gradient(pos_m, t, axis=0)
        a = np.gradient(v, t, axis=0)

        # POPRAWKA: tak jak w TIMDR-Radar-Module, ostatni gradient też
        # liczony względem rzeczywistego czasu t (oryginał mieszał
        # gradient "po czasie" z gradientem "po indeksie próbki").
        flow = np.gradient(v + a, t, axis=0)
        return flow

    # --- 2. Twist: nagłe zmiany kursu / wysokości ---
    def twist(self, track, angle_thresh=0.35, climb_rate_thresh_mps=15.0):
        """
        Wykrywa topologiczny twist:
          - nagłe zmiany kierunku (kursu) - angle_thresh w radianach
            (domyślnie 0.35 rad ~= 20 stopni)
          - nagłe zmiany prędkości pionowej - climb_rate_thresh_mps w m/s
            (domyślnie 15 m/s ~= 2950 ft/min; typowy komercyjny odrzutowiec
            wznosi się z prędkością rzędu 5-15 m/s, więc to próg "stromego"
            wznoszenia/zniżania - dostosuj do typu statku powietrznego,
            to nie jest zwalidowana wartość ATC)

        Zwraca słownik z indeksami punktów: direction_twist, altitude_twist.
        """
        tr = self._validate(track)
        lat, lon, alt, t = tr[:, 0], tr[:, 1], tr[:, 2], tr[:, 3]
        lat0 = float(np.mean(lat))

        xy = self._project_local_xy(lat, lon, lat0)
        v = np.gradient(xy, t, axis=0)

        # POPRAWKA (bug zawijania kąta, identyczny jak w TIMDR-Radar-Module):
        # kurs bliski 180°/-180° (lot na południe) powoduje, że arctan2
        # przeskakuje między wartościami blisko +pi i -pi. Naiwny
        # np.gradient() na takim ciągu kątów daje skoki ~360° zamiast
        # rzeczywistej zmiany kursu o kilka stopni. Zweryfikowano
        # empirycznie: lot na południe (kurs ~180°) z rzeczywistym
        # wahnięciem kursu ~0.1-0.6° dawał 4 fałszywe alarmy "twist" na
        # 5 punktów przy naiwnym gradiencie; 0 po poprawce. Poprawka:
        # np.unwrap() na ciągu kątów PRZED różniczkowaniem.
        angles_unwrapped = np.unwrap(np.arctan2(v[:, 1], v[:, 0]))
        dtheta = np.gradient(angles_unwrapped)

        # POPRAWKA (bug jednostek prędkości pionowej): oryginalny kod liczył
        # np.gradient(alt) bez `t`, czyli surową różnicę wysokości MIĘDZY
        # PRÓBKAMI, a nie prędkość pionową. Dla identycznej fizycznej
        # dynamiki lotu (te same 1000m różnicy wysokości) próg "50"
        # oznaczał co innego przy próbkowaniu co 10s i co 30s - dawał
        # dokładnie te same surowe liczby (200-300) niezależnie od
        # interwału próbkowania, mimo że rzeczywista prędkość pionowa
        # różniła się 3-krotnie (20-30 m/s vs 6.7-10 m/s). Poprawka:
        # np.gradient(alt, t) daje prawdziwą prędkość pionową w m/s,
        # niezależną od częstotliwości próbkowania.
        climb_rate = np.gradient(alt, t)

        twist_dir = np.where(np.abs(dtheta) > angle_thresh)[0]
        twist_alt = np.where(np.abs(climb_rate) > climb_rate_thresh_mps)[0]

        return {"direction_twist": twist_dir, "altitude_twist": twist_alt}

    # --- 3. TRM-reduction: stabilizacja toru lotu ---
    def trm_reduce(self, track):
        """
        TRM: prosta redukcja szumu toru lotu (średnia krocząca 3-punktowa
        na [lat, lon, alt]). Pierwszy i ostatni punkt pozostają bez zmian.

        Uwaga: uśrednianie surowych stopni lat/lon jest matematycznie
        bezpieczne (w przeciwieństwie do gradientu/kursu nie wymaga
        korekty cos(lat)), o ile tor NIE przecina antimeridianu -
        `_validate()` już to sprawdza.
        """
        tr = self._validate(track)
        pos = tr[:, :3]

        smooth = pos.copy()
        for i in range(1, len(pos) - 1):
            smooth[i] = (pos[i - 1] + pos[i] + pos[i + 1]) / 3.0

        return smooth

    # --- 4. Predykcja trajektorii ---
    def predict(self, track, steps=5):
        """
        Prosta predykcja pozycji [lat, lon, alt] na `steps` kroków w przód,
        zakładając lokalnie stałe przyspieszenie ekstrapolowane w lokalnym
        układzie metrycznym, a następnie rzutowane z powrotem na
        lat/lon. Ostrzeżenie: to ekstrapolacja czysto kinematyczna z
        ostatnich kilku próbek - nie modeluje planu lotu, wiatru ani
        intencji pilota/ATC. Wiarygodna tylko na krótkim horyzoncie i dla
        lotu bez gwałtownych manewrów.
        """
        tr = self._validate(track)
        lat, lon, alt, t = tr[:, 0], tr[:, 1], tr[:, 2], tr[:, 3]
        lat0 = float(np.mean(lat))

        xy = self._project_local_xy(lat, lon, lat0)
        pos_m = np.column_stack([xy, alt])

        v = np.gradient(pos_m, t, axis=0)
        a = np.gradient(v, t, axis=0)

        dt = t[-1] - t[-2]  # bezpieczne: _validate gwarantuje dt > 0

        pred_m = []
        p = pos_m[-1].copy()
        v0 = v[-1].copy()
        a0 = a[-1].copy()
        for _ in range(steps):
            v0 = v0 + a0 * dt
            p = p + v0 * dt
            pred_m.append(p.copy())
        pred_m = np.array(pred_m)

        lat_pred, lon_pred = self._inverse_local_xy(
            pred_m[:, :2], lat0, lat_ref=lat[0], lon_ref=lon[0]
        )
        return np.column_stack([lat_pred, lon_pred, pred_m[:, 2]])

    # --- diagnostyka pomocnicza (kurs / prędkość / prędkość pionowa) ---
    def diagnostics(self, track):
        """
        Zwraca per-punkt: kurs (stopnie, 0-360 od północy), prędkość
        względem ziemi (m/s i węzły), prędkość pionową (m/s i ft/min).
        Przydatne np. do wyświetlenia obok flag z twist().
        """
        tr = self._validate(track)
        lat, lon, alt, t = tr[:, 0], tr[:, 1], tr[:, 2], tr[:, 3]
        lat0 = float(np.mean(lat))

        xy = self._project_local_xy(lat, lon, lat0)
        v = np.gradient(xy, t, axis=0)
        climb_rate_mps = np.gradient(alt, t)

        # kurs mierzony od polnocy (0 st = polnoc, 90 st = wschod)
        course_deg = (np.rad2deg(np.arctan2(v[:, 0], v[:, 1]))) % 360.0
        ground_speed_mps = np.linalg.norm(v, axis=1)

        return {
            "course_deg": course_deg,
            "ground_speed_mps": ground_speed_mps,
            "ground_speed_kt": ground_speed_mps * 1.943844,
            "climb_rate_mps": climb_rate_mps,
            "climb_rate_fpm": climb_rate_mps * 196.850394,
        }
