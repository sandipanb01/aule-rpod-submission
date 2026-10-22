"""Time systems: UTC, TAI, TT (~TDB) and UT1, anchored at a UTC epoch.

Internally every time is ``t`` = SI seconds elapsed since the scenario UTC
epoch (float64, ~0.1 us resolution over a month). The constructor verifies
that no leap second occurs inside the analysis span, so ``t`` maps linearly
onto TAI and TT; ERFA two-part Julian dates are produced on demand.

Earth-orientation parameters (UT1-UTC, polar motion) come from the IERS
finals2000A file bundled with ``astropy-iers-data`` (measured values for
July 2026, i.e. not predictions).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone

import erfa
import numpy as np

from .constants import SECONDS_PER_DAY, TT_MINUS_TAI_S

ARCSEC = np.pi / (180.0 * 3600.0)


def _parse_iso(s: str) -> datetime:
    dt = datetime.fromisoformat(s.replace("Z", ""))
    return dt.replace(tzinfo=None)


@dataclass
class TimeSystem:
    epoch_utc: str
    span_end_utc: str
    eop: bool = True
    jd1: float = field(init=False)
    jd2: float = field(init=False)
    tai_minus_utc_s: float = field(init=False)

    def __post_init__(self) -> None:
        e = _parse_iso(self.epoch_utc)
        self._epoch_dt = e
        self.jd1, self.jd2 = erfa.dtf2d(
            "UTC", e.year, e.month, e.day, e.hour, e.minute, e.second + e.microsecond * 1e-6
        )
        self.tai_minus_utc_s = float(erfa.dat(e.year, e.month, e.day, 0.0))
        end = _parse_iso(self.span_end_utc)
        dat_end = float(erfa.dat(end.year, end.month, end.day, 0.0))
        if dat_end != self.tai_minus_utc_s:
            raise ValueError("A leap second occurs inside the analysis span; "
                             "piecewise TAI-UTC handling would be required.")
        self._span_s = (end - e).total_seconds()
        self._load_eop()

    # ------------------------------------------------------------------ EOP
    def _load_eop(self) -> None:
        days = np.arange(-2.0, self._span_s / SECONDS_PER_DAY + 3.0, 1.0)
        self._eop_t = days * SECONDS_PER_DAY
        if not self.eop:
            self._dut1 = np.zeros_like(days)
            self._xp = np.zeros_like(days)
            self._yp = np.zeros_like(days)
            return
        from astropy.utils import iers  # local import: only needed here
        iers.conf.auto_download = False
        tab = iers.IERS_A.open(iers.IERS_A_FILE)
        jd2 = self.jd2 + days
        dut1, st1 = tab.ut1_utc(self.jd1, jd2, return_status=True)
        xp, yp, st2 = tab.pm_xy(self.jd1, jd2, return_status=True)
        self.eop_status = {"ut1_utc": np.asarray(st1).tolist(), "pm": np.asarray(st2).tolist()}
        self._dut1 = np.asarray(dut1.to_value("s"), dtype=float)
        self._xp = np.asarray(xp.to_value("arcsec"), dtype=float) * ARCSEC
        self._yp = np.asarray(yp.to_value("arcsec"), dtype=float) * ARCSEC

    def dut1_s(self, t):
        return np.interp(t, self._eop_t, self._dut1)

    def polar_motion_rad(self, t):
        return np.interp(t, self._eop_t, self._xp), np.interp(t, self._eop_t, self._yp)

    # ------------------------------------------------------- Julian dates
    def jd_utc(self, t):
        return self.jd1, self.jd2 + np.asarray(t, dtype=float) / SECONDS_PER_DAY

    def jd_tt(self, t):
        t = np.asarray(t, dtype=float)
        return self.jd1, self.jd2 + (t + self.tai_minus_utc_s + TT_MINUS_TAI_S) / SECONDS_PER_DAY

    def jd_ut1(self, t):
        t = np.asarray(t, dtype=float)
        return self.jd1, self.jd2 + (t + self.dut1_s(t)) / SECONDS_PER_DAY

    # ------------------------------------------------------- conversions
    def t_from_iso(self, s: str) -> float:
        return (_parse_iso(s) - self._epoch_dt).total_seconds()

    def iso(self, t, decimals: int = 3) -> np.ndarray:
        """ISO-8601 UTC strings (valid because no leap second in span)."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        base = np.datetime64(self._epoch_dt.isoformat(), "us")
        us = np.round(t * 1e6).astype("int64").astype("timedelta64[us]")
        out = np.datetime_as_string(base + us, unit="us")
        if decimals < 6:
            cut = 20 + decimals if decimals > 0 else 19
            out = np.array([s[:cut] for s in out])
        return np.char.add(out, "Z")

    def datetime(self, t) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, dtype=float))
        base = np.datetime64(self._epoch_dt.isoformat(), "us")
        return base + np.round(t * 1e6).astype("int64").astype("timedelta64[us]")

    def ut1_hours_of_day(self, t: float) -> float:
        """UT1 time of day in hours (mean solar time at Greenwich)."""
        j1, j2 = self.jd_ut1(t)
        frac = (j1 + j2 - 0.5) % 1.0
        return float(frac * 24.0)
