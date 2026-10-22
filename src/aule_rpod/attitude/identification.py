"""Closed-form identification of the target attitude law.

Model (regular precession, expressed in the reference/LVLH frame):

    R_LVLH<-B(t) = Rot(a, w_p t) . R0 . Rot(b, w_s t)

a   : precession axis, fixed in LVLH          w_p : precession rate
b   : spin axis, fixed in the body            w_s : spin rate
R0  : attitude at the first sample            nutation = angle(a, R0 b)

This family contains every torque-free axisymmetric motion and every STK
'precessing spin' profile; an inertially fixed spin seen from LVLH is the
special case a = orbit normal, w_p = -n. A Levenberg-Marquardt fit to all
samples (geodesic residuals) recovers the parameters; a residual at machine
precision proves the law is exact, which (i) removes interpolation error
between the 60-s samples and (ii) explains the target's motion physically.

Initial guesses come from the data themselves: the per-step rotation axis
precesses about `a` in the reference frame at w_p and about `b` in the body
frame at -w_s, so their phase rates seed the fit (no manual tuning).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

from .stk_io import AttitudeProfile


def _unit(th, ph):
    return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])


def _sph(v):
    v = v / np.linalg.norm(v)
    return np.arccos(np.clip(v[2], -1, 1)), np.arctan2(v[1], v[0])


@dataclass
class PrecessionLaw:
    t0_s: float
    axis_ref: np.ndarray       # a (LVLH)
    rate_p_rad_s: float        # w_p
    axis_body: np.ndarray      # b (body)
    rate_s_rad_s: float        # w_s
    R0: Rotation
    rms_residual_rad: float = np.nan
    max_residual_rad: float = np.nan

    def rotation(self, t) -> Rotation:
        tt = np.atleast_1d(np.asarray(t, dtype=float)) - self.t0_s
        return (Rotation.from_rotvec(np.outer(self.rate_p_rad_s * tt, self.axis_ref))
                * self.R0
                * Rotation.from_rotvec(np.outer(self.rate_s_rad_s * tt, self.axis_body)))

    @property
    def nutation_deg(self) -> float:
        return float(np.degrees(np.arccos(np.clip(np.dot(self.axis_ref, self.R0.apply(self.axis_body)), -1, 1))))

    def summary(self) -> dict:
        return {
            "precession_axis_lvlh": np.round(self.axis_ref, 12).tolist(),
            "precession_rate_deg_s": float(np.degrees(self.rate_p_rad_s)),
            "precession_period_s": float(2 * np.pi / abs(self.rate_p_rad_s)),
            "spin_axis_body": np.round(self.axis_body, 12).tolist(),
            "spin_rate_deg_s": float(np.degrees(self.rate_s_rad_s)),
            "spin_period_s": float(2 * np.pi / abs(self.rate_s_rad_s)),
            "nutation_deg": self.nutation_deg,
            "R0_lvlh_from_body": self.R0.as_matrix().round(12).tolist(),
            "rms_residual_rad": float(self.rms_residual_rad),
            "max_residual_rad": float(self.max_residual_rad),
        }


def identify_precession_law(profile: AttitudeProfile) -> PrecessionLaw:
    t = profile.times_s - profile.times_s[0]
    R = Rotation.from_quat(profile.quat_xyzw)
    dt = np.median(np.diff(t))
    # --- data-driven initial guess ------------------------------------------
    w_ref = (R[1:] * R[:-1].inv()).as_rotvec()
    w_body = (R[:-1].inv() * R[1:]).as_rotvec()
    a0 = w_ref.mean(0)
    a0 /= np.linalg.norm(a0)
    b0 = w_body.mean(0)
    b0 /= np.linalg.norm(b0)

    def phase_rate(w, axis):
        e1 = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(axis, e1)
        ph = np.unwrap(np.arctan2(w @ e2, w @ e1))
        return np.polyfit(t[:-1], ph, 1)[0]

    wp0 = phase_rate(w_ref, a0)           # reference-frame wobble rate = +w_p
    ws0 = -phase_rate(w_body, b0)         # body-frame cone rate = -w_s

    def model(p, tt):
        a = _unit(p[0], p[1])
        b = _unit(p[2], p[3])
        return (Rotation.from_rotvec(np.outer(p[4] * tt, a)) * Rotation.from_rotvec(p[6:9])
                * Rotation.from_rotvec(np.outer(p[5] * tt, b)))

    def resid(p, RR, tt):
        return (model(p, tt).inv() * RR).as_rotvec().ravel()

    p0 = np.r_[_sph(a0), _sph(b0), wp0, ws0, R[0].as_rotvec()]
    n_seed = min(len(t), 600)
    fit = least_squares(resid, p0, args=(R[:n_seed], t[:n_seed]), method="lm",
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    fit = least_squares(resid, fit.x, args=(R, t), method="lm",
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=50000)
    err = np.linalg.norm(fit.fun.reshape(-1, 3), axis=1)
    p = fit.x
    return PrecessionLaw(
        t0_s=float(profile.times_s[0]), axis_ref=_unit(p[0], p[1]), rate_p_rad_s=float(p[4]),
        axis_body=_unit(p[2], p[3]), rate_s_rad_s=float(p[5]), R0=Rotation.from_rotvec(p[6:9]),
        rms_residual_rad=float(np.sqrt(np.mean(err ** 2))), max_residual_rad=float(err.max()))
