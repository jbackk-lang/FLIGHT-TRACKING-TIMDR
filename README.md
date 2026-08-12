# FLIGHT-TRACKING-TIMDR

Śledzenie trajektorii 3D dla modeli latających (drony, RC) — filtr Kalmana
(pozycja + prędkość w x/y/z) plus niezależny detektor krzywizny i torsji
trajektorii, oparty na pomyśle z repo `jbackk-lang/THE`, ale poprawiony po
rygorystycznym sprawdzeniu numerycznym (patrz niżej).

Samodzielny projekt — nie zależy od `RADAR-TRACKING-TIMDR` (2D, pojazdy
naziemne), architektura tylko wzorowana na tamtym repo.

---

## ⚠️ Status danych: SYNTETYCZNE

**Ten projekt na razie nie jest zwalidowany na prawdziwych logach lotu.**
`data/synthetic_flights.csv` (i skrypt `data/generate_synthetic_flights.py`,
który go generuje) to sztucznie wygenerowane trajektorie — linia prosta,
płaski zakręt, wznosząca się spirala, postój — z dodanym realistycznym, ale
też sztucznym szumem czujnika (sigma_xy=0.3m, sigma_z=0.5m, niezmierzone,
tylko orientacyjne rzędy wielkości).

Walidacja krzywizny/torsji jest za to mocniejsza niż w `RADAR-TRACKING-TIMDR`
w jednym sensie: dla idealnej matematycznej helisy krzywizna i torsja mają
**dokładny wzór analityczny** (patrz niżej), więc można sprawdzić zgodność
co do ułamka procenta, zamiast tylko korelacji z przybliżonym proxy jak przy
prawdziwych trasach GPS. Ale to nadal nie jest dowód, że model radzi sobie
z prawdziwym szumem czujnika lotu, manewrami pilota, wiatrem itd. Jeśli
podrzucisz prawdziwe logi (log lotu drona/RC z lat/lon/wysokość w czasie),
dorobię walidację na nich, tak jak przy `real_trips_sample.csv`.

---

## Dlaczego to jest osobny projekt, a nie rozszerzenie RADAR-TRACKING-TIMDR

Ten sam pomysł (krzywizna trajektorii jako sygnał manewru) już raz przeszedł
przez rygorystyczne odrzucenie i poprawkę w RADAR-TRACKING-TIMDR (2D, patrz
tamtejszy README, sekcja "Krzywizna trajektorii..."). Tam padły dwa wnioski:

1. Torsja nie ma sensu w 2D (płaszczyzna ruchu jest z definicji cała
   przestrzenią, więc "wychodzenie z płaszczyzny" jest zawsze zerowe).
2. Próba obejścia tego przez sztuczną "głębokość percepcyjną"
   `dz = f(dx, dy)` została odrzucona — deterministyczna funkcja
   istniejących danych nie wnosi nowej informacji, i dawała fałszywe
   sygnały na zwykłej linii prostej pod kątem do osi.

Modele latające to jedyny sensowny przypadek na realne rozszerzenie do 3D —
**pod warunkiem że Z pochodzi z prawdziwego czujnika** (barometr, GPS
altitude, radar wysokościowy), a nie z fabrykowanej funkcji x/y. Stąd osobne
repo: to jest właściwe miejsce na 3D, którego RADAR-TRACKING-TIMDR (naziemny,
bez czujnika wysokości) nigdy mieć nie będzie.

---

## Drugi błąd znaleziony po drodze: oryginalny wzór torsji "THE" nie zbiega do prawdy

Przy sprawdzaniu numerycznym na helisie o znanej analitycznie torsji, wzór
z oryginalnego pseudokodu THE:

```
tau = dot(cross(D_t-2, D_t-1), D_t) / G^2
```

(gdzie D to znormalizowane wektory kierunku) **systematycznie zbiega do
zera** zamiast do prawdziwej wartości torsji, gdy zagęszcza się próbkowanie
trajektorii. To osobny błąd matematyczny od problemu z fabrykowanym Z —
ujawnia się dopiero przy teście zbieżności z wzorem analitycznym, którego
nie było możliwe zrobić na prawdziwych, zaszumionych trasach GPS (nie ma
tam znanej "prawdziwej" torsji do porównania).

Poprawny wzór (użyty w `core/curvature_detector_3d.py`) liczy krzywiznę i
torsję ze standardowych wzorów różniczkowych na krzywą parametryczną,
estymując prędkość/przyspieszenie/szarpnięcie (v, a, j) skończonymi
różnicami z 4 kolejnych pozycji:

```
kappa = |v x a| / |v|^3
tau   = det(v, a, j) / |v x a|^2
```

### Walidacja na helisie o znanym wzorze analitycznym

Helisa: x=r·cos(t), y=r·sin(t), z=c·t ma dokładne wzory
`kappa = r/(r²+c²)`, `tau = c/(r²+c²)` (stałe wzdłuż całej krzywej).

| dt (gęstość próbkowania) | błąd kappa | błąd tau |
|---|---|---|
| 0.5 | 3.06% | 3.94% |
| 0.1 | 0.12% | 0.15% |
| 0.02 | 0.005% | 0.006% |

Błąd maleje z gęstszym próbkowaniem — to jest właściwa zbieżność, w
przeciwieństwie do oryginalnego wzoru, który przy zagęszczaniu zbiegał
do zera niezależnie od prawdziwej wartości.

Dodatkowe kontrole: płaski okrąg (bez wznoszenia) daje poprawnie
`kappa=1/r` oraz `tau=0` (poprawnie wykrywa płaskość); linia prosta w
dowolnym kierunku 3D daje `kappa=tau=0`; wynik jest niezmienniczy względem
obrotu układu współrzędnych wokół dowolnej osi 3D (nie tylko wokół Z).

---

## Ten sam problem co zawsze: dzielenie blisko zera przy niskiej prędkości

Krzywizna dzieli przez `|v|^3`, torsja przez `|v x a|^2` — przy postoju lub
bardzo wolnym locie szum czujnika (kilkadziesiąt cm) jest wzmacniany do
fałszywych pików, dokładnie jak w JRegulatorze, Helix-Astro T1 i
2D-wersji tego detektora. Naprawa: próg `min_speed`, poniżej którego wynik
jest bramkowany (`gated=True`, wartości 0.0) zamiast liczony.

W wersji 2D próg (3.0 m/krok) dobrano empirycznie z realnych tras GPS. Tu,
bo działamy na danych syntetycznych, domyślna wartość (`1.0`) jest punktem
startowym do skalibrowania na realnych logach lotu, nie zwalidowaną stałą.
Test regresyjny (`test_ungated_hover_would_blow_up_regression_guard`)
dowodzi, że bez progu ten sam szum czujnika daje krzywiznę >1.0 zamiast
poprawnych ~0 — żeby nikt przypadkiem nie usunął zabezpieczenia myśląc,
że jest zbędne.

### Drugi, osobny przypadek tego samego wzorca błędu: torsja na niemal prostej trajektorii

Znaleziony przy okazji sprawdzania wariantu pseudokodu z repo THE
(`THE_GEO_PRO_4D_Radar`), który zabezpieczał dzielenie w torsji warunkiem
`if cross_norm == 0`. Ten kod miał tę samą lukę co pierwotna wersja tu:
`cross_norm < 1e-12` chroni tylko przed dosłownym zerem, nie przed
wzmocnieniem szumu, gdy trajektoria jest niemal (ale nie dokładnie)
prostoliniowa — np. `v=(1,0,0), a=(1,1e-6,0)` (typowy szum kierunku na
"prostym" odcinku lotu) daje `cross_norm=1e-6`, nie łapie się w próg
`1e-12`, a `tau = dot(cross_va,j)/cross_norm**2` wychodzi rzędu 1e6
zamiast ~0.

Naprawa: bramkowanie torsji na podstawie samej krzywizny `kappa` (już
obliczonej, fizycznie sensownej wielkości), nie na `cross_norm` wprost —
`if kappa < min_curvature: tau = 0`. Domyślne `min_curvature=1e-4`, tak
jak `min_speed`, to punkt startowy wymagający kalibracji na realnych
danych. Test regresyjny:
`test_near_straight_noisy_direction_gives_zero_torsion_not_blowup`.

---

## Struktura

```
core/
  kalman_filter_3d.py       # filtr Kalmana stałej prędkości, stan 6D (x,y,z,vx,vy,vz)
  curvature_detector_3d.py  # krzywizna + torsja z progiem min_speed
  flight_tracker_3d.py      # łączy oba w jeden interfejs .update(x,y,z)
data/
  generate_synthetic_flights.py  # generator danych testowych (patrz wyżej)
  synthetic_flights.csv          # wygenerowane dane, 4 segmenty x ~100 kroków
tests/
  test_kalman_filter_3d.py
  test_curvature_detector_3d.py   # w tym test zgodności z helisą analityczną
  test_flight_tracker_3d.py       # integracja na syntetycznych danych
```

## Użycie

```python
from core.flight_tracker_3d import FlightTracker3D

tracker = FlightTracker3D(min_curvature_speed=1.0, min_curvature=1e-4)
for x, y, z in strumien_pomiarow:
    (est_x, est_y, est_z), krzywizna = tracker.update(x, y, z)
    if not krzywizna.gated and krzywizna.torsion > prog:
        ...  # podejrzana zmiana płaszczyzny lotu (np. korkociąg, spirala)
```

## Testy

```
python3 -m pytest tests/ -v
```

17 testów bazowego pipeline'u (Kalman + krzywizna/torsja), opisanych
wyżej. Poniższe cztery moduły dokładają kolejnych 26 (43 razem).

---

## Klasyfikator manewrów (`core/maneuver_classifier.py`)

Klasyfikuje typ manewru z krzywizny/torsji: prosto / zakręt / spiralizacja
/ korkociąg / gwałtowny manewr / nieznany (gated). Progi kalibrowane
empirycznie na 6 oznaczonych segmentach syntetycznych (`data/synthetic_flights.csv`
+ nowy `data/synthetic_maneuvers.csv`, dodaje `tight_corkscrew` i
`sharp_maneuver`), na wygładzonej medianą-z-5-próbek krzywiźnie -- pojedyncza
próbka jest zbyt zaszumiona do progowania wprost (patrz docstring modułu).

Rozróżnienie "gwałtowny manewr" vs trwały zakręt/spiralę opiera się na
CZASIE TRWANIA podwyższonej krzywizny, nie jej wielkości -- klasyfikator
online nie może wiedzieć z góry, czy początek jest krótkim epizodem czy
trwałym manewrem, więc każdy oznacza najpierw jako "gwałtowny manewr" i
dopiero po `confirm_samples` (domyślnie 8) próbkach "awansuje" do
zakrętu/spirali/korkociągu.

Walidacja (majority-vote na każdym z 6 segmentów, `tests/test_maneuver_classifier.py`):
straight_climb 78% prosto, level_turn 55% zakręt, ascending_helix 90%
spiralizacja, tight_corkscrew 81% korkociąg, hover 100% nieznany,
sharp_maneuver 83% prosto (+ epizody gwałtownego manewru dokładnie w
oknie faktycznego skrętu). level_turn ma najsłabszy wynik (55%, reszta
myli się jako spiralizacja) -- uczciwie udokumentowane w kodzie, nie
ukryte: promień 30 daje kappa blisko progu szumu przy tym poziomie
symulowanego szumu GPS/barometru.

## Predyktor trajektorii (`core/trajectory_predictor.py`)

Przewiduje przyszłą pozycję całkując równania Freneta-Serreta (RK4)
zamiast zwykłej ekstrapolacji liniowej -- zakłada, że bieżąca
krzywizna/torsja/prędkość utrzymują się na horyzoncie predykcji.
Zwalidowane na analitycznej helisie: błąd < 1e-4 przy delta_t do 2s (50
podkroków RK4/predykcję); na tej samej helisie bije ekstrapolację
liniową o ponad 5x na dalszym horyzoncie (`tests/test_trajectory_predictor.py`).
Dla toru lokalnie prostoliniowego (kappa~0) poprawnie spada do zwykłej
ekstrapolacji liniowej.

## Adaptacyjny Kalman (`core/adaptive_kalman_filter_3d.py`)

Rozszerza `KalmanFilter3D` o: (1) adaptacyjne Q przez innovation-based
adaptive estimation (IAE, klasyczna technika Mehra) -- Q rośnie
`q_boost_factor`-krotnie, gdy średnia energia ostatnich innowacji
przekracza próg, samowystarczalne, nie wymaga zewnętrznego sygnału
krzywizny; (2) opcjonalny `measurement_quality` (0,1] podawany przy
`update()`, skalujący R -- domyślnie 1.0 (brak informacji), więc bez
zewnętrznego sygnału jakości zachowuje się identycznie jak statyczne R
(świadomie "martwe" bez prawdziwego źródła SNR, ale przetestowane z
jawnie wstrzykiwaną syntetyczną jakością).

Walidacja (`tests/test_adaptive_kalman_filter_3d.py`): na segmencie
sharp_maneuver adaptacyjne Q daje 48% mniejszy błąd śledzenia w oknie
manewru niż statyczne Q, przy porównywalnym (lekko lepszym) błędzie na
odcinkach prostych. Wstrzyknięta niska `measurement_quality` podczas
sztucznego burstu szumu redukuje błąd o ~24% względem ignorowania jej.

## Wizualizator (`tools/visualizer.py`)

Generuje PNG na segment: trajektoria 3D (pomiar surowy / referencja
prawdziwa / estymata Kalmana), krzywizna/torsja/helikalność w czasie z
zaznaczonym bramkowaniem, pasek klasyfikacji manewru w czasie, błąd
śledzenia vs prawda. Uruchomienie: `python3 tools/visualizer.py` ->
`tools/output/<segment>.png`.

Ciekawa, uczciwie pozostawiona obserwacja z wizualizacji ascending_helix:
błąd estymaty Kalmana bywa WIĘKSZY niż błąd surowego pomiaru -- model
stałej prędkości jest systematycznie niedopasowany do zakrzywionego toru
(lag), to nie błąd w kodzie, tylko realne ograniczenie modelu stałej
prędkości na krzywiznie (stąd sens adaptacyjnego Kalmana i predyktora
opartego o krzywiznę wyżej).

## Fuzja radar+IMU (`core/imu_radar_fusion_3d.py`)

Luźno sprzężona fuzja: przyspieszenie z IMU jako wejście sterujące w
predykcji Kalmana (zamiast modelu stałej prędkości), pozycja z radaru
koryguje w update. Wymagało nowych syntetycznych danych IMU
(`data/generate_synthetic_imu.py` -> `data/synthetic_imu.csv`:
przyspieszenie prawdziwe z różnicy centralnej x_true/y_true/z_true +
szum biały 0.15 m/s² + dryft/bias random walk -- jawnie oznaczone jako
syntetyczne).

Kluczowa poprawka znaleziona przy budowie (opisana w pełni w docstringu
modułu): Q NIE może być dowolną małą stałą -- musi być wyprowadzone z
faktycznego szumu IMU propagowanego przez macierz sterowania
(`Q = B @ Sigma_accel @ B^T`), inaczej filtr ufa własnej zaszumionej
predykcji tak samo jak modelowi stałej prędkości i fuzja wypada GORZEJ
niż sam radar.

Walidacja z UCZCIWIE dopasowaną `measurement_var=0.3` dla obu wariantów
(zgodną z faktycznym szumem symulacji pozycji, nie domyślnym 5.0 z
KalmanFilter3D, które jest zawyżone względem tego konkretnego modelu
szumu): w oknie manewru fuzja daje błąd śr. 0.688 vs 0.870 radar-only
(-21%), ale na odcinku prostym 0.507 vs 0.357 (+42%, GORZEJ) --
udokumentowany, realny kompromis luźnej fuzji INS/GPS przy tej jakości
IMU, nie ukryta wada. Martwa nawigacja z samego IMU (bez radaru,
`ImuDeadReckoning3D`) dryfuje bez ograniczeń (błąd rośnie >5x od
początku do końca segmentu w teście) -- uzasadnienie, dlaczego fuzja (a
nie sam IMU) jest potrzebna.

---

## Testy (pełny pakiet)

```
python3 -m pytest tests/ -v
```

43 testy, wszystkie przechodzą (17 bazowych: 2 Kalman + 11
krzywizna/torsja + 4 integracja, + 8 klasyfikator manewrów + 5
predyktor trajektorii + 6 adaptacyjny Kalman + 2 wizualizator + 5 fuzja
radar+IMU).

## Licencja

MIT — patrz `LICENSE`.
