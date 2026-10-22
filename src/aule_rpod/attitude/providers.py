"""Attitude providers: continuous-time R_LVLH<-Body(t) from the 60-s samples.

slerp : spherical linear interpolation (constant body rate between samples) --
        the baseline; faithful to the delivered data, C0-continuous.
nlerp : component-wise linear interpolation + normalisation, i.e. STK's default
        'Lagrange, order 1' attitude interpolation.
law   : the identified closed-form precession law (exact between samples).
All three agree to <0.22 deg in docking-port direction (see validation report).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation, Slerp

from .identification import identify_precession_law
from .stk_io import AttitudeProfile


class SlerpAttitude:
    name = "SLERP of 60-s samples"

    def __init__(self, profile: AttitudeProfile):
        self.profile = profile
        self._slerp = Slerp(profile.times_s, Rotation.from_quat(profile.quat_xyzw))

    def rotation(self, t) -> Rotation:
        t = np.clip(np.atleast_1d(np.asarray(t, dtype=float)), self.profile.times_s[0], self.profile.times_s[-1])
        return self._slerp(t)


class NlerpAttitude:
    name = "NLERP (STK Lagrange order 1)"

    def __init__(self, profile: AttitudeProfile):
        self.profile = profile

    def rotation(self, t) -> Rotation:
        tq, q = self.profile.times_s, self.profile.quat_xyzw
        t = np.clip(np.atleast_1d(np.asarray(t, dtype=float)), tq[0], tq[-1])
        k = np.clip(np.searchsorted(tq, t, side="right") - 1, 0, len(tq) - 2)
        f = ((t - tq[k]) / (tq[k + 1] - tq[k]))[:, None]
        qi = (1 - f) * q[k] + f * q[k + 1]
        return Rotation.from_quat(qi / np.linalg.norm(qi, axis=1, keepdims=True))


class LawAttitude:
    name = "identified closed-form law"

    def __init__(self, profile: AttitudeProfile, law=None):
        self.profile = profile
        self.law = law or identify_precession_law(profile)

    def rotation(self, t) -> Rotation:
        return self.law.rotation(t)


def make_attitude_provider(kind: str, profile: AttitudeProfile):
    return {"slerp": SlerpAttitude, "nlerp": NlerpAttitude, "law": LawAttitude}[kind](profile)
