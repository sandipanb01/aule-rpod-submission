"""Analytic first-order secular J2 propagator.

Elements a, e, i are constant; Omega, omega and M drift linearly
(Vallado 2013 Eqs. 9-38/9-39/9-41; Brouwer 1959 first-order secular terms):

    n      = sqrt(mu/a^3)
    k      = J2 (Re/p)^2
    dOm/dt = -3/2 n k cos i
    dw/dt  =  3/4 n k (5 cos^2 i - 1)
    dM/dt  =  n [1 + 3/4 k sqrt(1-e^2) (3 cos^2 i - 1)]

The position/velocity at any time follow from the two-body relations of the
drifted elements. This is the formulation of STK's 'J2 Perturbation'
propagator ("secular drift only ... based upon Keplerian mean elements") with
its default 'Osculating' ellipse option, i.e. the entered a, e, i are used
directly. Closed form => exact, fast, and smooth for event root-finding.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .elements import kep2cart, true_from_mean


@dataclass(frozen=True)
class SecularJ2Orbit:
    a_km: float
    e: float
    i_rad: float
    raan0_rad: float
    argp0_rad: float
    M0_rad: float
    mu: float
    re_km: float
    j2: float
    name: str = "secular J2 (analytic)"

    @property
    def rates(self):
        n = np.sqrt(self.mu / self.a_km ** 3)
        p = self.a_km * (1.0 - self.e ** 2)
        k = self.j2 * (self.re_km / p) ** 2
        ci = np.cos(self.i_rad)
        omdot = -1.5 * n * k * ci
        wdot = 0.75 * n * k * (5.0 * ci ** 2 - 1.0)
        mdot = n * (1.0 + 0.75 * k * np.sqrt(1.0 - self.e ** 2) * (3.0 * ci ** 2 - 1.0))
        return n, omdot, wdot, mdot

    @property
    def nodal_period_s(self) -> float:
        _, _, wdot, mdot = self.rates
        return 2.0 * np.pi / (wdot + mdot)

    def elements(self, t):
        t = np.asarray(t, dtype=float)
        _, omdot, wdot, mdot = self.rates
        raan = self.raan0_rad + omdot * t
        argp = self.argp0_rad + wdot * t
        M = self.M0_rad + mdot * t
        return raan, argp, M

    def state_teoe(self, t):
        raan, argp, M = self.elements(t)
        nu = true_from_mean(M, self.e)
        return kep2cart(self.a_km, self.e, self.i_rad, raan, argp, nu, self.mu)

    def arg_latitude(self, t):
        raan, argp, M = self.elements(t)
        return np.mod(argp + true_from_mean(M, self.e), 2 * np.pi)
