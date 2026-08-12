import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from core.kalman_filter_3d import KalmanFilter3D


def test_tracks_constant_velocity_3d_target():
    random.seed(1)
    kf = KalmanFilter3D(process_var=0.01, measurement_var=2.0)

    true_pos = [0.0, 0.0, 0.0]
    velocity = (3.0, -1.0, 0.5)

    errors = []
    for _ in range(80):
        true_pos = [true_pos[i] + velocity[i] for i in range(3)]
        measured = tuple(true_pos[i] + random.gauss(0, 1.4) for i in range(3))
        est = kf.update(measured)
        err = sum((est[i] - true_pos[i]) ** 2 for i in range(3)) ** 0.5
        errors.append(err)

    # after the initial transient, estimate should track closely despite noise
    assert sum(errors[-20:]) / 20 < 3.0

    vx, vy, vz = kf.velocity_estimate
    assert abs(vx - velocity[0]) < 0.5
    assert abs(vy - velocity[1]) < 0.5
    assert abs(vz - velocity[2]) < 0.5


def test_predict_only_advances_without_measurement():
    kf = KalmanFilter3D(initial_position=(0.0, 0.0, 0.0))
    kf.update((1.0, 1.0, 1.0))
    kf.update((2.0, 2.0, 2.0))
    before = kf.x.copy()
    est = kf.predict()
    assert est != (before[0, 0], before[1, 0], before[2, 0])  # state advanced
