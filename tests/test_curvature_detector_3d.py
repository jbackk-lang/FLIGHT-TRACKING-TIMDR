import sys, os, math, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.curvature_detector_3d import CurvatureDetector3D


def test_first_three_updates_are_gated_no_history():
    det = CurvatureDetector3D(min_speed=0.0)
    for _ in range(3):
        r = det.update(0.0, 0.0, 0.0)
        assert r.gated is True
        assert r.curvature == 0.0
        assert r.torsion == 0.0


def test_straight_line_any_3d_direction_gives_zero_kappa_and_tau():
    det = CurvatureDetector3D(min_speed=0.1)
    direction = (2.0, -3.0, 1.5)
    pts = [(direction[0] * t, direction[1] * t, direction[2] * t) for t in range(8)]
    results = [det.update(x, y, z) for x, y, z in pts]
    for r in results[3:]:
        assert abs(r.curvature) < 1e-9
        assert abs(r.torsion) < 1e-9
        assert r.gated is False


def test_flat_circle_gives_correct_curvature_and_zero_torsion():
    """A planar circle (constant altitude) has an exact known curvature
    1/r and MUST have zero torsion -- torsion measures departure from a
    plane, and this path never leaves its plane."""
    det = CurvatureDetector3D(min_speed=0.1)
    r = 10.0
    dt = 0.05
    pts = [(r * math.cos(t * dt), r * math.sin(t * dt), 50.0) for t in range(30)]
    results = [det.update(x, y, z) for x, y, z in pts]
    mid = results[15]
    assert not mid.gated
    assert abs(mid.curvature - 1.0 / r) / (1.0 / r) < 0.01  # within 1%
    assert abs(mid.torsion) < 1e-6


def test_ascending_helix_matches_analytical_curvature_and_torsion():
    """Exact ground truth exists for a helix: kappa = r/(r^2+c^2),
    tau = c/(r^2+c^2). This is the strongest possible validation -- not
    correlation with a noisy proxy, but agreement with closed-form math."""
    det = CurvatureDetector3D(min_speed=0.1)
    r, c = 5.0, 1.0
    dt = 0.02
    n = 400
    pts = [(r * math.cos(t * dt), r * math.sin(t * dt), c * t * dt) for t in range(n)]
    results = [det.update(x, y, z) for x, y, z in pts]

    kappa_true = r / (r**2 + c**2)
    tau_true = c / (r**2 + c**2)

    mid = results[n // 2]
    assert not mid.gated
    assert abs(mid.curvature - kappa_true) / kappa_true < 0.01
    assert abs(mid.torsion - tau_true) / tau_true < 0.02


def test_hover_with_sensor_noise_is_gated():
    random.seed(0)
    det = CurvatureDetector3D(min_speed=3.0)  # comfortably above the ~0.3-0.5m sensor noise floor
    results = []
    for _ in range(20):
        x = random.gauss(0, 0.3)
        y = random.gauss(0, 0.3)
        z = 40.0 + random.gauss(0, 0.5)
        results.append(det.update(x, y, z))
    for r in results[3:]:
        assert r.gated is True
        assert r.curvature == 0.0
        assert r.torsion == 0.0


def test_ungated_hover_would_blow_up_regression_guard():
    """Same rationale as the 2D detector: proves *why* the gate exists.
    Without it, sensor noise at near-zero speed is amplified into a huge
    spurious curvature (division by speed^3)."""
    random.seed(0)
    det = CurvatureDetector3D(min_speed=0.0)
    results = []
    for _ in range(20):
        x = random.gauss(0, 0.3)
        y = random.gauss(0, 0.3)
        z = 40.0 + random.gauss(0, 0.5)
        results.append(det.update(x, y, z))
    max_kappa = max(r.curvature for r in results[3:])
    assert max_kappa > 1.0


def _rotation_matrix(axis, theta):
    norm = sum(a * a for a in axis) ** 0.5
    axis = tuple(a / norm for a in axis)
    ax, ay, az = axis
    c, s = math.cos(theta), math.sin(theta)
    C = 1 - c
    return [
        [c + ax*ax*C,      ax*ay*C - az*s, ax*az*C + ay*s],
        [ay*ax*C + az*s,   c + ay*ay*C,    ay*az*C - ax*s],
        [az*ax*C - ay*s,   az*ay*C + ax*s, c + az*az*C],
    ]


def _apply(mat, p):
    return tuple(sum(mat[i][j] * p[j] for j in range(3)) for i in range(3))


def test_rotation_invariance_arbitrary_3d_axis():
    """Curvature and torsion of a physical path must not depend on the
    orientation of the sensor's coordinate frame. Rotate the same helix
    around an arbitrary (non axis-aligned) 3D axis and confirm the
    curvature/torsion sequence is unchanged."""
    r, c = 4.0, 0.8
    dt = 0.05
    path = [(r * math.cos(t * dt), r * math.sin(t * dt), c * t * dt) for t in range(40)]

    def run(points):
        det = CurvatureDetector3D(min_speed=0.01)
        return [det.update(*p) for p in points]

    base = run(path)

    mat = _rotation_matrix((0.4, 0.7, 0.3), 1.1)
    rotated_path = [_apply(mat, p) for p in path]
    rotated = run(rotated_path)

    for a, b in zip(base[3:], rotated[3:]):
        assert abs(a.curvature - b.curvature) < 1e-6
        assert abs(a.torsion - b.torsion) < 1e-6


def test_near_straight_noisy_direction_gives_zero_torsion_not_blowup():
    """Regresja: znaleziony przy okazji sprawdzania pseudokodu
    THE_GEO_PRO_4D_Radar blad -- stary warunek `cross_norm < 1e-12`
    chronil tylko przed dzieleniem przez doslowne zero, nie przed
    wzmocnieniem szumu gdy trajektoria jest niemal (ale nie dokladnie)
    prosta. Ten test wymusza male, ale nie zerowe odchylenie kierunku
    (typowy szum sensora) i sprawdza, ze torsja zostaje wyzerowana
    przez bramkowanie na kappa, zamiast eksplodowac do rzedu 1e6."""
    det = CurvatureDetector3D(min_speed=0.0)
    det.update(0.0, 0.0, 0.0)
    det.update(1.0, 0.0, 0.0)
    det.update(2.0, 0.0, 0.0)
    r = det.update(3.0, 1e-6, 0.0)  # male odchylenie zamiast idealnej prostej
    assert r.gated is False
    assert abs(r.torsion) < 1e-6


def test_min_curvature_gate_does_not_affect_real_helix():
    """Upewnij sie, ze nowe bramkowanie na kappa nie psuje prawdziwej,
    dobrze zakrzywionej trajektorii (kappa >> domyslny prog)."""
    det = CurvatureDetector3D(min_speed=0.0)
    r, c = 5.0, 1.0
    kappa_true = r / (r * r + c * c)
    tau_true = c / (r * r + c * c)
    dt = 0.02
    for k in range(4):
        t = 2.0 + k * dt
        res = det.update(r * math.cos(t), r * math.sin(t), c * t)
    assert res.gated is False
    assert abs(res.curvature - kappa_true) / kappa_true < 0.01
    assert abs(res.torsion - tau_true) / tau_true < 0.01


def test_min_curvature_rejects_negative():
    try:
        CurvatureDetector3D(min_curvature=-0.1)
        assert False, "powinno rzucic ValueError"
    except ValueError:
        pass


def test_min_speed_rejects_negative():
    try:
        CurvatureDetector3D(min_speed=-1.0)
        assert False, "should have raised"
    except ValueError:
        pass
