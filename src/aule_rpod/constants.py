"""Physical constants with provenance.

Earth constants are read from ``config/mission.yaml`` at run time (so a
sensitivity study can change them); the values below are the defaults and the
references they come from.
"""

# --- Earth (WGS-84 / EGM96) ---------------------------------------------------
MU_EARTH_KM3_S2 = 398600.4418          # EGM96 GM [km^3/s^2]
R_EARTH_KM = 6378.137                  # WGS-84 equatorial radius [km]
F_EARTH = 1.0 / 298.257223563          # WGS-84 flattening [-]
J2_EARTH = 1.08262668e-3               # EGM96, -C20 (unnormalised) [-]

# --- Sun ----------------------------------------------------------------------
AU_KM = 149_597_870.700                # IAU 2012 Resolution B2 [km]
R_SUN_KM = 695_700.0                   # IAU 2015 Resolution B3 nominal radius [km]

# --- Time ---------------------------------------------------------------------
SECONDS_PER_DAY = 86_400.0
TT_MINUS_TAI_S = 32.184                # exact, by definition
TROPICAL_YEAR_DAYS = 365.2421897       # mean tropical year (J2000) [days]

DEG = 3.141592653589793 / 180.0
