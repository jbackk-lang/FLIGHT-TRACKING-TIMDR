"""
demo_conflict_map.py -- wizualizacja conflict_alert() na mapie (3 gotowe
przyklady), uruchamiana z menu w run.bat.

Rysuje dla kazdego przykladu: historyczny tor obu samolotow (linia ciagla),
przewidywana trajektorie (linia przerywana, ta sama kinematyka co
predict()), i jesli conflict_alert() zglasza konflikt - czerwony znacznik
w przewidywanym punkcie najwiekszego zblizenia, z podpisem separacji.

To ilustracja dzialania funkcji, NIE nowe obliczenia bezpieczenstwa -
liczby na wykresie pochodza wprost z conflict_alert()/predict(), ktore sa
juz przetestowane w test_conflict_alert.py.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")  # zapis do pliku, bez wymogu okna/ekranu
import matplotlib.pyplot as plt

from timdr_flight import TIMDRFlight

timdr = TIMDRFlight()


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


def _rysuj_przyklad(ax, nazwa, track_a, track_b, lookahead=180.0, pred_steps=20):
    wynik = timdr.conflict_alert(track_a, track_b, lookahead_seconds=lookahead)

    pred_a = timdr.predict(track_a, steps=pred_steps)
    pred_b = timdr.predict(track_b, steps=pred_steps)

    ax.plot(track_a[:, 1], track_a[:, 0], "o-", color="tab:blue", label="Samolot A (historia)")
    ax.plot(track_b[:, 1], track_b[:, 0], "o-", color="tab:orange", label="Samolot B (historia)")
    ax.plot(pred_a[:, 1], pred_a[:, 0], "--", color="tab:blue", alpha=0.5, label="A (predykcja)")
    ax.plot(pred_b[:, 1], pred_b[:, 0], "--", color="tab:orange", alpha=0.5, label="B (predykcja)")

    if wynik["conflict"]:
        dt_track_a = track_a[1, 3] - track_a[0, 3]
        dt_track_b = track_b[1, 3] - track_b[0, 3]
        i_a = min(int(round(wynik["time_to_conflict_s"] / dt_track_a)), pred_steps - 1)
        i_b = min(int(round(wynik["time_to_conflict_s"] / dt_track_b)), pred_steps - 1)
        ax.plot(pred_a[i_a, 1], pred_a[i_a, 0], "rx", markersize=14, markeredgewidth=3)
        ax.plot(pred_b[i_b, 1], pred_b[i_b, 0], "rx", markersize=14, markeredgewidth=3)
        status = (f"KONFLIKT za {wynik['time_to_conflict_s']:.0f}s "
                  f"(poziomo {wynik['min_horizontal_nm']:.2f} NM, "
                  f"pionowo {wynik['min_vertical_ft']:.0f} ft)")
        kolor_tytulu = "tab:red"
    else:
        status = (f"Brak konfliktu (min. separacja: "
                  f"{wynik['min_horizontal_nm']:.2f} NM poziomo / "
                  f"{wynik['min_vertical_ft']:.0f} ft pionowo)")
        kolor_tytulu = "tab:green"

    ax.set_title(f"{nazwa}\n{status}", fontsize=10, color=kolor_tytulu)
    ax.set_xlabel("długość geogr.")
    ax.set_ylabel("szerokość geogr.")
    ax.legend(fontsize=7, loc="best")
    ax.grid(alpha=0.3)


def main():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    # Przyklad 1: kurs kolizyjny, ta sama wysokosc -> KONFLIKT
    a1 = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b1 = _prosty_tor(50.0, 20.5, 10000.0, heading_deg=270, speed_mps=200)
    _rysuj_przyklad(axes[0], "Przykład 1: kurs kolizyjny, ta sama wysokość", a1, b1)

    # Przyklad 2: kurs kolizyjny poziomo, ale 5000ft roznicy wysokosci -> bezpiecznie
    a2 = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b2 = _prosty_tor(50.0, 20.5, 11524.0, heading_deg=270, speed_mps=200)
    _rysuj_przyklad(axes[1], "Przykład 2: separacja pionowa 5000 ft", a2, b2)

    # Przyklad 3: lot rownolegly, 10 NM odstepu -> bezpiecznie
    a3 = _prosty_tor(50.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    b3 = _prosty_tor(50.0 + 18520.0 / 111320.0, 19.5, 10000.0, heading_deg=90, speed_mps=200)
    _rysuj_przyklad(axes[2], "Przykład 3: lot równoległy, 10 NM odstępu", a3, b3)

    fig.suptitle("FLIGHT-TRACKING-TIMDR: conflict_alert() - 3 przykłady", fontsize=13)
    fig.tight_layout()

    out_path = "conflict_alert_przyklady.png"
    fig.savefig(out_path, dpi=130)
    print(f"Zapisano mapę z przykładami: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
