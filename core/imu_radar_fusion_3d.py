"""
imu_radar_fusion_3d.py
------------------------
Luzno sprzezona (loosely-coupled) fuzja radar+IMU: przyspieszenie z IMU
jest wejsciem sterujacym (control input) w kroku predykcji filtra
Kalmana (zamiast zwyklego modelu stalej predkosci), a pozycja z radaru
koryguje predykcje w kroku aktualizacji -- standardowa architektura
INS/GPS, tu nazwana zgodnie z terminologia repo "radar" (patrz README).

Dlaczego to ma sens (i kiedy):
  - Radar/GPS: dokladna dlugoterminowo (brak dryfu), ale rzadsza/
    zaszumiona pozycja, i -- co wazniejsze przy manewrach -- model
    stalej predkosci w zwyklym KalmanFilter3D NIE WIE o przyspieszeniu
    do momentu, az zobaczy je posrednio w kolejnych pomiarach pozycji
    (opoznienie/lag, patrz test_adaptive_kalman -- to samo zjawisko).
  - IMU: mierzy przyspieszenie BEZPOSREDNIO i natychmiast, wiec model
    ruchu "wie" o manewrze od razu -- ale samo calkowanie przyspieszenia
    (bez korekty pozycja) dryfuje bez ograniczen (bias, szum sumuje sie
    kwadratowo w pozycji przez podwojne calkowanie -- patrz
    test_imu_only_drifts_without_correction).
  - Fuzja: przyspieszenie z IMU daje predykcje "swiadoma manewru" bez
    czekania na kolejne pomiary radaru (mniejszy lag niz zwykly Kalman
    podczas manewru), a radar okresowo koryguje pozycje, zapobiegajac
    nieograniczonemu dryfowi samego IMU.

Model stanu: [x,y,z,vx,vy,vz], dt=1 (jak reszta repo). Predykcja z
przyspieszeniem jako control input:
    x_pred  = x + vx*dt + 0.5*ax*dt^2
    vx_pred = vx + ax*dt
(analogicznie y,z). Korekta pozycja z radaru -- identyczna jak w
KalmanFilter3D (H mierzy tylko pozycje).

Kalibracja Q (WAZNE, znaleziona empirycznie przy budowie tego modulu):
pierwsza wersja uzywala stalej, malej Q=0.01*I (jak w KalmanFilter3D) i
fuzja wypadala GORZEJ niz sam radar -- nawet podczas manewru. Powod:
filtr nie "wiedzial", ze jego wlasna predykcja (napedzana zaszumionym
IMU) jest niepewna, wiec przewazal ja nad korekta z radaru. Poprawka:
Q = B @ Sigma_accel @ B^T, czyli szum procesu jest WYPROWADZONY z
faktycznego szumu przyspieszenia (accel_noise_std) propagowanego przez
te sama macierz sterowania B, ktora go wstrzykuje do stanu -- filtr
poprawnie "wie", jak bardzo nie ufac wlasnej predykcji.

Walidacja (patrz tests/test_imu_radar_fusion_3d.py), UCZCIWE porownanie
z TA SAMA measurement_var=0.3 dla obu wariantow (dopasowana do
faktycznego szumu symulacji pozycji: sigma_xy=0.3 -> war.=0.09,
sigma_z=0.5 -> war.=0.25 z generate_synthetic_flights.py -- wczesniejszy
domyslny measurement_var=5.0 z KalmanFilter3D byl znacznie zawyzony
wzgledem tego konkretnego modelu szumu):

    segment sharp_maneuver, okno manewru (probki 38-54):
        radar-only Kalman:   blad śr. 0.870
        fuzja radar+IMU:     blad śr. 0.688   (-21%)

    ten sam segment, odcinek prosty (probki 0-34):
        radar-only Kalman:   blad śr. 0.357
        fuzja radar+IMU:     blad śr. 0.507   (+42%, GORZEJ)

To jest UCZCIWY, typowy dla luznej fuzji INS/GPS kompromis, nie ukryta
wada: IMU o tej jakosci szumu (accel_noise_std=0.15 m/s^2) pomaga
wyraznie podczas manewru (natychmiastowa informacja o przyspieszeniu,
mniejszy lag niz model stalej predkosci), ale na spokojnym odcinku
dokladany szum przyspieszeniomierza (podwojnie calkowany do pozycji)
jest wiekszym zrodlem bledu niz sam szum pozycji radaru. W realnym
systemie rozwiazuje sie to adaptacyjnie (ufaj IMU wiecej tylko podczas
wykrytego manewru -- patrz core/adaptive_kalman_filter_3d.py, ktory
robi dokladnie to, tylko bez IMU jako wejscia) -- polaczenie obu technik
to naturalny nastepny krok, nie zrobiony tutaj (zakres tego modulu to
podstawowa, uczciwie zwalidowana fuzja luzna).

ImuDeadReckoning3D (martwa nawigacja z samego IMU, bez radaru) jest
dolaczona jako punkt odniesienia pokazujacy, DLACZEGO fuzja (a nie sam
IMU) jest potrzebna -- bez korekty pozycja jej blad rosnie w
przyblizeniu kwadratowo z czasem (patrz
test_imu_only_drifts_without_correction: rzedu dziesiatek metrow po
kilkunastu krokach, setek po calym segmencie).
"""
import numpy as np


class ImuRadarFusion3D:
    def __init__(
        self,
        dt: float = 1.0,
        accel_noise_std: float = 0.15,
        process_var: float = 0.001,
        measurement_var: float = 0.3,
        initial_position=(0.0, 0.0, 0.0),
        initial_velocity=(0.0, 0.0, 0.0),
        initial_uncertainty: float = 100.0,
    ):
        if dt <= 0:
            raise ValueError("dt musi byc > 0")
        self.dt = dt
        x0, y0, z0 = initial_position
        vx0, vy0, vz0 = initial_velocity
        self.x = np.array([[x0], [y0], [z0], [vx0], [vy0], [vz0]])
        self.P = np.eye(6) * initial_uncertainty

        self.F = np.eye(6)
        self.F[0, 3] = dt
        self.F[1, 4] = dt
        self.F[2, 5] = dt

        # macierz sterowania: przyspieszenie (ax,ay,az) -> zmiana [x,y,z,vx,vy,vz]
        self.B = np.zeros((6, 3))
        self.B[0, 0] = 0.5 * dt ** 2
        self.B[1, 1] = 0.5 * dt ** 2
        self.B[2, 2] = 0.5 * dt ** 2
        self.B[3, 0] = dt
        self.B[4, 1] = dt
        self.B[5, 2] = dt

        self.H = np.zeros((3, 6))
        self.H[0, 0] = 1.0
        self.H[1, 1] = 1.0
        self.H[2, 2] = 1.0

        self.R = np.eye(3) * measurement_var

        # KLUCZOWE dla poprawnej fuzji: Q NIE jest dowolna mala stala, tylko
        # jest WYPROWADZONE z szumu przyspieszenia propagowanego przez B
        # (Q = B * Sigma_accel * B^T). Bez tego filtr "ufa" wlasnej
        # predykcji napedzanej zaszumionym IMU dokladnie tak samo jak
        # predykcji ze stalej predkosci -- i systematycznie przewaza
        # zaszumiony sygnal IMU nad korekta z radaru (sprawdzone
        # empirycznie: ze stala process_var=0.01 fuzja wypadala GORZEJ niz
        # sam radar, bo filtr nie "wiedzial", ze jego wlasna predykcja jest
        # niepewna). process_var to dodatkowy, mniejszy margines na
        # niemodelowana dynamike (bias IMU, blad dt=1 itp.).
        sigma_accel_sq = accel_noise_std ** 2
        Sigma_accel = np.eye(3) * sigma_accel_sq
        Q_imu = self.B @ Sigma_accel @ self.B.T
        self.Q = Q_imu + np.eye(6) * process_var

    def predict(self, accel):
        ax, ay, az = accel
        u = np.array([[ax], [ay], [az]])
        self.x = self.F @ self.x + self.B @ u
        self.P = self.F @ self.P @ self.F.T + self.Q
        return (float(self.x[0, 0]), float(self.x[1, 0]), float(self.x[2, 0]))

    def correct(self, radar_position):
        zx, zy, zz = radar_position
        z_vec = np.array([[zx], [zy], [zz]])
        y = z_vec - self.H @ self.x
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)
        self.x = self.x + K @ y
        self.P = (np.eye(6) - K @ self.H) @ self.P
        return (float(self.x[0, 0]), float(self.x[1, 0]), float(self.x[2, 0]))

    def update(self, accel, radar_position):
        """Pelny krok: predykcja z IMU, potem korekta z radaru."""
        self.predict(accel)
        return self.correct(radar_position)

    @property
    def velocity_estimate(self):
        return (float(self.x[3, 0]), float(self.x[4, 0]), float(self.x[5, 0]))


class ImuDeadReckoning3D:
    """Martwa nawigacja z samego IMU, BEZ korekty pozycja -- referencyjny
    "zly przypadek", zeby pokazac, dlaczego fuzja jest potrzebna
    (patrz docstring modulu i test_imu_only_drifts_without_correction).
    """

    def __init__(self, dt: float = 1.0, initial_position=(0.0, 0.0, 0.0), initial_velocity=(0.0, 0.0, 0.0)):
        self.dt = dt
        self.position = list(initial_position)
        self.velocity = list(initial_velocity)

    def update(self, accel):
        ax, ay, az = accel
        for i, a in enumerate((ax, ay, az)):
            self.position[i] += self.velocity[i] * self.dt + 0.5 * a * self.dt ** 2
            self.velocity[i] += a * self.dt
        return tuple(self.position)
