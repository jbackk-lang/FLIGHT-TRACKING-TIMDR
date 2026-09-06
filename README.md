# THE
TIMDR Hyperflow Engine (THE)

🧠 THE — rdzeń kodowy (pseudokod)

To jest minimalna, koncepcyjna struktura THE. Sekcje 1-7 poniżej to
**język/metafora**, nie zwalidowany numerycznie kod — nie ma tu wielkości
fizycznych do sprawdzenia względem wartości analitycznej, więc nie były
testowane w ten sposób. Sekcja "THE-GEO PRO" niżej jest inna: opisuje
konkretne wielkości geometryczne (krzywizna, skręt trajektorii), które
*da się* sprawdzić względem wartości analitycznych — i po sprawdzeniu
okazało się, że oryginalny wzór na skręt był błędny. Historia tej
poprawki jest opisana w sekcji "THE-GEO PRO → THE-GEO PRO 4D" niżej.

Nie jest to implementacja w żadnym języku — to jest język THE.

## 1. Strumień (S-Layer)
Strumień nie przechowuje wartości — tylko zmianę.

```
class Strumien:
    zmiana = 0

    def update(dane):
        zmiana = dane - zmiana

    def gradient():
        return zmiana

    def pulse():
        return abs(zmiana)
```
To jest odpowiednik „procesora”, ale bez CPU.

## 2. Topologia (T-Layer)
Topologia nie jest grafem danych — jest grafem przepływu.

```
class Topologia:
    wezly = []
    krawedzie = []

    def route(strumien):
        return krawedzie[strumien]

    def reshape():
        reorganizuj_polaczenia(wezly, krawedzie)
```
To jest odpowiednik „kolejek”, ale bez kolejek.

## 3. Przepływ (F-Layer)
Przepływ nie jest schedulerem — jest kierunkiem zmiany.

```
class Przeplyw:
    def direction(zmiana):
        if zmiana > 0: return "UP"
        if zmiana < 0: return "DOWN"
        return "STABLE"

    def velocity(zmiana):
        return abs(zmiana)

    def reorganize():
        dostosuj_kierunki()
```
To jest odpowiednik „schedulerów”, ale bez schedulerów.

## 4. Stabilność (C-Layer)
Stabilność nie jest kontrolą błędów — jest filtracją percepcji.

```
class Stabilnosc:
    def signal(zmiana):
        return abs(zmiana) < prog

    def noise(zmiana):
        return abs(zmiana) > prog

    def coherence(strumien):
        return strumien.gradient() < limit
```
To jest odpowiednik „kontroli błędów”, ale bez błędów.

## 🔥 5. THE Hyperflow Loop — główna pętla percepcyjna
Zamiast CPU → scheduler → proces → wątek → blokada → kolejka, masz:
strumień → topologia → przepływ → stabilność

```
while True:
    S.update(dane)
    T.reshape()
    F.reorganize()
    C.coherence(S)
```
To jest pętla percepcyjna, nie obliczeniowa.

## 🧬 6. THE — przepływ helikalny (opcjonalny moduł)

```
class HelikalnyPrzeplyw:
    def cycle(zmiana):
        return sin(zmiana)

    def linear(zmiana):
        return zmiana

    def combine():
        return cycle(zmiana) + linear(zmiana)
```
To jest pipeline bez pipeline.

## 🚀 7. THE — minimalny system

```
S = Strumien()
T = Topologia()
F = Przeplyw()
C = Stabilnosc()

while True:
    S.update(input)
    T.reshape()
    F.reorganize()
    C.coherence(S)
```
To jest pełny THE w 12 liniach pseudokodu.

🧠 **Najprostsza definicja THE kodu**: cztery klasy (strumień, topologia,
przepływ, stabilność) połączone w pętlę percepcyjną.

---

## THE-GEO PRO → THE-GEO PRO 4D: historia walidacji

Poniższe sekcje dotyczą jednej konkretnej, mierzalnej wielkości:
geometrii trajektorii punktu w przestrzeni (kierunek, krzywizna, skręt).
W przeciwieństwie do sekcji 1-7 wyżej, te wzory dają się sprawdzić —
istnieje krzywa analityczna (helisa) o znanej, dokładnej krzywiźnie i
skręcie, więc każdą proponowaną formułę można porównać z prawdziwą
wartością.

### Oryginalny wzór (poniżej, sekcja "THE-GEO PRO") — sprawdzony i odrzucony

Oryginalna `Torsion.compute` liczyła skręt z różnic **kierunków
jednostkowych** `D_t2, D_t1, D_t` (patrz kod niżej). Sprawdzone na
czystej, bezszumowej helisie analitycznej `x=r·cos(t), y=r·sin(t), z=c·t`
(r=5, c=1, prawdziwe τ=0.038462):

| dt | τ (stary wzór) | błąd |
|---|---|---|
| 0.5 | 0.003507 | 90.9% |
| 0.1 | 0.000724 | 98.1% |
| 0.02 | 0.000145 | 99.6% |
| 0.01 | 0.000073 | 99.8% |
| 0.005 | 0.000036 | 99.9% |

Błąd **rośnie** w stronę 100% zamiast maleć do zera przy zagęszczaniu
próbkowania — to znaczy, że wzór nie zbiega do prawdziwego skrętu, tylko
do zera. Wzór jest błędny niezależnie od tego, jak gęste są dane.

### Poprawiony wzór — THE-GEO PRO 4D

Zamiast różnic kierunków jednostkowych, poprawna wersja liczy krzywiznę
i skręt z pochodnych **prędkości / przyspieszenia / szarpnięcia**
(v, a, j — różnice skończone z 4 kolejnych punktów, gdzie t = parametr
czasowy krzywej, stąd "4D" = x,y,z,t):

```
v  = p_t  - p_t1
v1 = p_t1 - p_t2
a  = v - v1
a1 = v1 - (p_t2 - p_t3)
j  = a - a1

speed = |v|
if speed < min_speed:
    return {gated: True, curvature: 0, torsion: 0, helical: 0}

kappa = |v x a| / |v|^3
tau   = det(v, a, j) / |v x a|^2
H     = sqrt(kappa^2 + tau^2)
```

Sprawdzone na tej samej helisie (r=5, c=1, κ=0.192308, τ=0.038462):

| dt | κ | błąd κ | τ | błąd τ |
|---|---|---|---|---|
| 0.5 | 0.186416 | 3.06% | 0.039977 | 3.94% |
| 0.1 | 0.192070 | 0.12% | 0.038521 | 0.15% |
| 0.02 | 0.192298 | 0.0049% | 0.038464 | 0.0062% |
| 0.01 | 0.192305 | 0.0012% | 0.038462 | 0.0015% |
| 0.005 | 0.192307 | 0.0003% | 0.038462 | 0.0004% |

Błąd maleje monotonicznie do ~0 — wzór faktycznie zbiega do prawdziwej
krzywizny i skrętu.

**Implementacja**: [`the_geo_pro_4d.py`](the_geo_pro_4d.py) —
prawdziwy, uruchamialny kod Python (nie pseudokod).

### Druga poprawka: bramkowanie oparte na cross_norm==0 nie wystarcza

Wersja pseudokodu nadesłana później pod nazwą `THE_GEO_PRO_4D_Radar`
zabezpieczała dzielenie przez `if cross_norm == 0`. To chroni tylko
przed dosłownym zerem/NaN — **nie** chroni przed wzmocnieniem szumu, gdy
trajektoria jest niemal (ale nie dokładnie) prosta:

```
v = (1, 0, 0);  a = (1, 1e-6, 0)   # typowy szum kierunku na "prostym" odcinku
cross_norm = 1e-6   # != 0, więc stary warunek przepuszcza dalej
tau = dot(cross_va, j) / cross_norm**2  →  rzędu 1e6 zamiast ~0
```

Ten sam wzorzec błędu co przy `min_speed`/`min_step_m` (dzielenie przez
małą wartość wzmacnia szum), tylko na innym mianowniku. Poprawka:
bramkowanie torsji na podstawie krzywizny `kappa` (już obliczonej,
fizycznie sensownej wielkości), nie na `cross_norm` wprost —
`if kappa < min_curvature: tau = 0`. Domyślne `min_curvature=1e-4` to
punkt startowy, nie zwalidowana stała — jak `min_step_m` czy
`min_speed`, wymaga kalibracji na realnych zaszumionych danych 3D, gdy
się pojawią.

**Testy**: [`tests/test_the_geo_pro_4d.py`](tests/test_the_geo_pro_4d.py)
— 8 testów: zbieżność do wartości analitycznej, brak dzielenia przez
zero na linii prostej, bramkowanie `min_speed`, regresja na wzmacnianie
szumu przy niemal-prostej trajektorii (opisana wyżej), niezmienniczość
na obrót 3D wokół dowolnej osi (wzór Rodriguesa). Wszystkie przechodzą:

```
$ python3 -m unittest discover -s tests -v
...
Ran 10 tests in 0.003s
OK
```

### Trzecia poprawka: `THE_GEO_PRO_4D_Radar.py` było drugą, niezależną kopią tego samego wzoru

Znalezione przy pełnej ponownej inspekcji repo: `THE_GEO_PRO_4D_Radar.py`
miało WŁASNĄ implementację — treściowo identyczną z `the_geo_pro_4d.py`
(ta sama poprawka bramkowania na `kappa` z sekcji wyżej), ale jako
osobna kopia kodu, nie import ze wspólnego źródła. To dokładnie ten
mechanizm, który już raz spowodował błąd opisany w "Druga poprawka" wyżej
— dwie kopie tej samej logiki, jedna naprawiona niezależnie od drugiej.
Poprawka: `THE_GEO_PRO_4D_Radar.py` jest teraz cienkim re-eksportem
(`from the_geo_pro_4d import THE_GEO_PRO_4D as THE_GEO_PRO_4D_Radar`),
nie osobną implementacją — nie da się już, żeby te dwie nazwy cicho się
rozjechały. Test regresyjny na tożsamość obiektu funkcji (nie tylko
"daje ten sam wynik teraz"):
[`tests/test_the_geo_pro_4d_radar.py`](tests/test_the_geo_pro_4d_radar.py).

### Czwarta poprawka: wzmacnianie szumu przy potrójnym różniczkowaniu na realnych/krótkich sekwencjach

Znalezione przy audycie ekosystemu TIMDR pod kątem powtarzalnego wzorca
błędu ("Pattern A"): `THE_GEO_PRO_4D()` liczy `v/a/j` z **surowych,
kolejnych** punktów. Na gładkiej, gęsto próbkowanej krzywej analitycznej
(jak w testach zbieżności wyżej) to zbiega poprawnie. Ale repo siostrzane
`FLIGHT-TRACKING-TIMDR`, używające matematycznie identycznego wzoru,
zmierzyło bezpośrednio na zaszumionych danych GPS: surowe różnicowanie
wzmacnia szum pomiaru 10-40× (`std(tau)≈2.6` przy prawdziwym `tau≈0.065`).
Niezależnie, `GIA-TIMDR` potwierdziło ten sam mechanizm na krótkich
(~25-punktowych) realnych szeregach czasowych — krzywizna/torsja z
surowych różnic stają się nierozróżnialne od szumu.

**Poprawka**: nowa funkcja `THE_GEO_PRO_4D_sequence(points, smooth=True)`
— dla sekwencji dłuższych niż 4 punkty liczy `v/a/j` **bezpośrednio z
pochodnych lokalnie dopasowanego wielomianu** (prawdziwy filtr
różniczkujący Savitzky-Golay, bez zależności od scipy — własna
eliminacja Gaussa dla równań normalnych najmniejszych kwadratów), nie
"wygładź pozycje, potem policz proste różnice sąsiadów". **Ta pierwsza,
prostsza wersja została wypróbowana i odrzucona** — dawała WIĘKSZY błąd
niż brak wygładzania w ogóle, bo odejmowanie dwóch niezależnie
dopasowanych, zachodzących na siebie okien nie tłumi szumu wyższych
pochodnych tak, jak wzięcie pochodnej analitycznie z jednego dopasowania
na punkt. Zmierzone na zaszumionej helisie: średni błąd `tau` spada
zauważalnie (test `test_smoothed_sequence_has_lower_tau_error_than_raw`)
względem surowego różnicowania. `THE_GEO_PRO_4D()` (4-punktowa) pozostaje
niezmieniona — `smooth=False` w nowej funkcji odtwarza jej stare
zachowanie dokładnie, dla wstecznej kompatybilności.

**Testy**: [`tests/test_noise_amplification_fix.py`](tests/test_noise_amplification_fix.py)
— 4 nowe testy (redukcja błędu na zaszumionej helisie, brak regresji na
czystych danych, dokładna zgodność `smooth=False` ze starym
zachowaniem, czytelny błąd zamiast `IndexError` przy `polyorder<3`).

### Zastosowania (już zbudowane i zwalidowane w osobnych repo)

- **[RADAR-TRACKING-TIMDR](https://github.com/jbackk-lang/RADAR-TRACKING-TIMDR)**
  — wariant 2D (bez skrętu, tylko krzywizna) sprawdzony na prawdziwych
  danych GPS z 4 przejazdów: korelacja z realnymi manewrami 0.47-0.76
  (przy bramkowaniu `min_step_m=3.0`) vs -0.08..-0.00 bez bramkowania.
- **[FLIGHT-TRACKING-TIMDR](https://github.com/jbackk-lang/FLIGHT-TRACKING-TIMDR)**
  — pełny wariant 3D (ten opisany wyżej) użyty do śledzenia lotu na
  danych syntetycznych; `the_geo_pro_4d.py` w tym repo jest matematycznie
  identyczny z `core/curvature_detector_3d.py` w FLIGHT-TRACKING-TIMDR
  (ponownie zweryfikowane numerycznie przy tej inspekcji — 200 losowych
  zestawów punktów, max różnica |Δkappa|≈4e-17, |Δtau|≈3e-16, czyli
  szum zaokrągleń zmiennoprzecinkowych, nie realna różnica), różni się
  tylko interfejsem (dict zamiast dataclass, dodatkowo liczy `H`).

---

## THE-GEO PRO (oryginalny pseudokod, zachowany dla kontekstu)

Poniższe klasy to oryginalny pseudokod tego repo. `Torsion.compute` jest
**błędny** — patrz sekcja wyżej. Zachowany tu bez zmian dla
przejrzystości historii, nie do użycia.

### 1. Delta geometryczna (pełna różniczka)
```
class DeltaGeo:
    def compute(p_t, p_t1):
        dx = p_t.x - p_t1.x
        dy = p_t.y - p_t1.y
        dz = p_t.z - p_t1.z
        return (dx, dy, dz)
```

### 2. Gradient + kierunek
```
class GradientDir:
    def gradient(dx, dy, dz):
        return sqrt(dx*dx + dy*dy + dz*dz)

    def direction(dx, dy, dz, G):
        if G == 0: return (0,0,0)
        return (dx/G, dy/G, dz/G)
```

### 3. Krzywizna (curvature)
```
class Curvature:
    def compute(D_t, D_t1, G):
        diff = norm(D_t - D_t1)
        return diff / G
```

### 4. Skręt (torsion) — ⚠️ BŁĘDNY, patrz walidacja wyżej
```
class Torsion:
    def compute(D_t2, D_t1, D_t, G):
        cross_vec = cross(D_t2, D_t1)
        return dot(cross_vec, D_t) / (G*G)
```

### 5. Helikalność (spiralność)
```
class Helical:
    def compute(kappa, tau):
        return sqrt(kappa*kappa + tau*tau)
```

### 6. Przepływ geometryczny PRO
```
class FlowGeo:
    def compute(D, G, kappa, tau):
        return {
            "dir": D,
            "vel": G,
            "curv": kappa,
            "tors": tau
        }
```

### 7. Stabilność geometryczna PRO
```
class StabilityGeo:
    def dir_stab(D_t, D_t1):
        return dot(D_t, D_t1)

    def curv_stab(kappa):
        return 1 / (1 + kappa)

    def tors_stab(tau):
        return 1 / (1 + abs(tau))

    def helix_stab(H):
        return 1 / (1 + H)

    def total(Dstab, Ck, Ct, Ch):
        return Dstab * Ck * Ct * Ch
```

### 8. THE-GEO PRO — główna pętla percepcyjna (oryginalna, z błędnym skrętem)
```
def THE_GEO_PRO(p_t, p_t1, p_t2):
    dx, dy, dz = DeltaGeo.compute(p_t, p_t1)
    G = GradientDir.gradient(dx, dy, dz)
    D_t = GradientDir.direction(dx, dy, dz, G)

    dx1, dy1, dz1 = DeltaGeo.compute(p_t1, p_t2)
    G1 = GradientDir.gradient(dx1, dy1, dz1)
    D_t1 = GradientDir.direction(dx1, dy1, dz1, G1)

    kappa = Curvature.compute(D_t, D_t1, G)

    dx2, dy2, dz2 = DeltaGeo.compute(p_t2, p_t1)
    G2 = GradientDir.gradient(dx2, dy2, dz2)
    D_t2 = GradientDir.direction(dx2, dy2, dz2, G2)

    tau = Torsion.compute(D_t2, D_t1, D_t, G)   # <- błędny wzór, patrz wyżej

    H = Helical.compute(kappa, tau)
    F = FlowGeo.compute(D_t, G, kappa, tau)

    Dstab = StabilityGeo.dir_stab(D_t, D_t1)
    Ck = StabilityGeo.curv_stab(kappa)
    Ct = StabilityGeo.tors_stab(tau)
    Ch = StabilityGeo.helix_stab(H)
    C = StabilityGeo.total(Dstab, Ck, Ct, Ch)

    return {
        "delta": (dx, dy, dz),
        "gradient": G,
        "direction": D_t,
        "curvature": kappa,
        "torsion": tau,
        "helical": H,
        "flow": F,
        "stability": C
    }
```

**Użyj zamiast tego [`the_geo_pro_4d.py`](the_geo_pro_4d.py)** — poprawny,
przetestowany, zwalidowany na helisie analitycznej.
