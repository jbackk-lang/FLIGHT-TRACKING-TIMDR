"""
FlightTracker3D -- łączy KalmanFilter3D (wygładzanie pozycji) z opcjonalnym
CurvatureDetector3D (diagnostyka manewru: krzywizna + torsja) w jeden
prosty interfejs .update(x, y, z).

Samodzielny odpowiednik RadarTrackerCustom z RADAR-TRACKING-TIMDR, ale bez
zależności między repozytoriami -- ten projekt nie importuje niczego z
tamtego, żeby dało się go rozpakować i użyć osobno.

CurvatureDetector3D tu NIE reguluje filtra (nie ma tu odpowiednika
JRegulatora) -- liczy diagnostykę równolegle, na wygładzonych przez Kalmana
pozycjach (co dodatkowo tłumi szum pomiaru zanim trafi do liczenia
krzywizny/torsji).
"""
from core.kalman_filter_3d import KalmanFilter3D
from core.curvature_detector_3d import CurvatureDetector3D, CurvatureResult3D


class FlightTracker3D:
    def __init__(
        self,
        min_curvature_speed: float = 1.0,
        min_curvature: float = 1e-4,
        **kalman_kwargs,
    ):
        self.filter = KalmanFilter3D(**kalman_kwargs)
        self.curvature = CurvatureDetector3D(
            min_speed=min_curvature_speed, min_curvature=min_curvature
        )

    def update(self, x: float, y: float, z: float):
        est = self.filter.update((x, y, z))
        curv = self.curvature.update(*est)
        return est, curv
