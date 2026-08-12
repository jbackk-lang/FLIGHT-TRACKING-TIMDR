"""
3D constant-velocity Kalman filter, single (x, y, z) measurement per step.

Stan: [x, y, z, vx, vy, vz]^T (6-wymiarowy). Rozszerzenie 1D wersji z
RADAR-TRACKING-TIMDR (core/kalman_filter_custom.py) na trzy niezależne osie
-- model ruchu jest blokowo-diagonalny (każda oś x/y/z porusza się
niezależnie ze stałą prędkością + szum procesu), więc to naprawdę jest po
prostu trzy równoległe filtry 1D złożone w jeden stan, nie coś bardziej
egzotycznego.

Ta sama zastrzeżenie co w wersji 1D: F/Q są zbudowane dla stałego dt=1 na
wywołanie update(). Przy nieregularnych odstępach czasowych między
pomiarami trzeba przebudować F/Q z prawdziwym dt.
"""
import numpy as np


class KalmanFilter3D:
    def __init__(self, process_var: float = 0.01, measurement_var: float = 5.0,
                 initial_position=(0.0, 0.0, 0.0), initial_uncertainty: float = 100.0):
        x0, y0, z0 = initial_position
        self.x = np.array([[x0], [y0], [z0], [0.0], [0.0], [0.0]])
        self.P = np.eye(6) * initial_uncertainty

        # blokowo-diagonalny model stalej predkosci: kazda os (x,y,z) niezalezna
        self.F = np.eye(6)
        self.F[0, 3] = 1.0  # x += vx
        self.F[1, 4] = 1.0  # y += vy
        self.F[2, 5] = 1.0  # z += vz

        self.H = np.zeros((3, 6))
        self.H[0, 0] = 1.0
        self.H[1, 1] = 1.0
        self.H[2, 2] = 1.0

        self.R = np.eye(3) * measurement_var
        self.Q = np.eye(6) * process_var

    def update(self, z):
        zx, zy, zz = z
        z_vec = np.array([[zx], [zy], [zz]])

        # predykcja
        self.x = self.F @ self.x
        self.P = self.F @ self.P @ self.F.T + self.Q

        # korekcja
        y = z_vec - self.H @ self.x
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)

        self.x = self.x + K @ y
        self.P = (np.eye(6) - K @ self.H) @ self.P

        return (float(self.x[0, 0]), float(self.x[1, 0]), float(self.x[2, 0]))

    def predict(self):
        """Sama predykcja bez korekcji pomiarem -- coasting na modelu ruchu."""
        self.x = self.F @ self.x
        self.P = self.F @ self.P @ self.F.T + self.Q
        return (float(self.x[0, 0]), float(self.x[1, 0]), float(self.x[2, 0]))

    @property
    def velocity_estimate(self):
        return (float(self.x[3, 0]), float(self.x[4, 0]), float(self.x[5, 0]))
