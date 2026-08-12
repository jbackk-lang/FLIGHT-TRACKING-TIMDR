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

tracker = FlightTracker3D(min_curvature_speed=1.0)
for x, y, z in strumien_pomiarow:
    (est_x, est_y, est_z), krzywizna = tracker.update(x, y, z)
    if not krzywizna.gated and krzywizna.torsion > prog:
        ...  # podejrzana zmiana płaszczyzny lotu (np. korkociąg, spirala)
```

## Testy

```
python3 -m pytest tests/ -v
```

14 testów, wszystkie przechodzą: filtr Kalmana 3D (śledzenie stałej
prędkości), krzywizna/torsja (linia prosta w 3D, płaski okrąg, helisa vs
wzór analityczny, bramkowanie przy postoju + test regresyjny dowodzący
konieczności bramkowania, niezmienniczość względem obrotu wokół dowolnej
osi 3D), integracja na syntetycznych danych (filtr redukuje szum, spirala
daje wyższą krzywiznę/torsję niż linia prosta, postój jest bramkowany).

## Licencja

MIT — patrz `LICENSE`.
