import numpy as np

class TIMDRFlight:
    def __init__(self):
        pass

    def timdr_flow(self, track):
        tr = np.array(track)
        pos = tr[:, :3]
        t = tr[:, 3]
        v = np.gradient(pos, t, axis=0)
        a = np.gradient(v, t, axis=0)
        flow = np.gradient(v + a, axis=0)
        return flow

    def twist(self, track, angle_thresh=0.35, alt_thresh=50.0):
        tr = np.array(track)
        pos = tr[:, :3]
        latlon = pos[:, :2]
        alt = pos[:, 2]
        v = np.gradient(latlon, axis=0)
        angles = np.arctan2(v[:,1], v[:,0])
        dtheta = np.gradient(angles)
        dalt = np.gradient(alt)
        twist_dir = np.where(np.abs(dtheta) > angle_thresh)[0]
        twist_alt = np.where(np.abs(dalt) > alt_thresh)[0]
        return {"direction_twist": twist_dir, "altitude_twist": twist_alt}

    def trm_reduce(self, track):
        tr = np.array(track)
        pos = tr[:, :3]
        smooth = pos.copy()
        for i in range(1, len(pos)-1):
            smooth[i] = (pos[i-1] + pos[i] + pos[i+1]) / 3.0
        return smooth

    def predict(self, track, steps=5):
        tr = np.array(track)
        pos = tr[:, :3]
        t = tr[:, 3]
        v = np.gradient(pos, t, axis=0)
        a = np.gradient(v, t, axis=0)
        dt = t[-1] - t[-2]
        p = pos[-1]
        v0 = v[-1]
        a0 = a[-1]
        pred = []
        for _ in range(steps):
            v0 = v0 + a0 * dt
            p = p + v0 * dt
            pred.append(p.copy())
        return np.array(pred)
