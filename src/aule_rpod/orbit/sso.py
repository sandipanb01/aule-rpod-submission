"""Sun-synchronous orbit design.

Inclination: first-order secular J2 nodal rate equal to the mean motion of the
mean Sun (one revolution per tropical year):

    dOmega/dt = -3/2 n J2 (Re/p)^2 cos i = 2 pi / T_trop     (Vallado 2013, Eq. 9-41)

RAAN from LTAN: the local MEAN solar time at terrestrial longitude lambda is
UT1 + lambda/15 (UT1 is Greenwich mean solar time), so at the orbit epoch the
ascending node sits at  lambda_N = 15 deg/h * (LTAN - UT1_h).  Its right
ascension in the TEOE frame (CIO-based) is  Omega = ERA(epoch) + lambda_N.
Polar motion (~0.3") is neglected in this one-off node placement.
"""
from __future__ import annotations

import numpy as np

from ..constants import SECONDS_PER_DAY


def sso_inclination_rad(a_km: float, e: float, mu: float, re_km: float, j2: float,
                        tropical_year_days: float) -> float:
    n = np.sqrt(mu / a_km ** 3)
    p = a_km * (1.0 - e * e)
    omdot = 2.0 * np.pi / (tropical_year_days * SECONDS_PER_DAY)
    cos_i = -omdot / (1.5 * n * j2 * (re_km / p) ** 2)
    if abs(cos_i) > 1:
        raise ValueError("No sun-synchronous inclination exists for this semi-major axis")
    return float(np.arccos(cos_i))


def raan_from_ltan_teoe(ltan_hours: float, ut1_hours_of_day: float, era_rad: float) -> float:
    lam_node = np.radians(15.0 * (ltan_hours - ut1_hours_of_day))
    return float(np.mod(era_rad + lam_node, 2.0 * np.pi))
