"""Independent cross-validation with Orekit 13 (CS GROUP / CNES flight-dynamics library).

Orekit is used as a black-box reference for:
  1. Numerical J2 propagation (J2OnlyPerturbation + DormandPrince853) from the
     same osculating state as our Cowell model, in the same frozen-CIRF frame.
  2. GCRF -> ITRF (IERS 2010, full EOP) and topocentric elevation for every
     station (OneAxisEllipsoid WGS-84 + TopocentricFrame).
  3. Ground-station AOS/LOS events (ElevationDetector) vs our event finder.
  4. Earth-shadow events (EclipseDetector, penumbra-entry convention) using a
     Sun position provider backed by ERFA (Orekit's own JPL loader needs DE
     files that are not shipped here).
Orekit data: a minimal data directory built from the IERS finals2000A file and
the IERS/ERFA leap-second table (data/orekit-data).
"""
from __future__ import annotations

import numpy as np

_VM = None


def _init(data_dir: str):
    global _VM
    if _VM is None:
        import orekit_jpype as orekit
        _VM = orekit.initVM()
        from orekit_jpype.pyhelpers import setup_orekit_data
        setup_orekit_data(filenames=str(data_dir), from_pip_library=False)


def run(sc, data_dir: str):
    _init(data_dir)
    from jpype import JImplements, JOverride
    from org.hipparchus.geometry.euclidean.threed import Vector3D
    from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
    from org.orekit.bodies import GeodeticPoint, OneAxisEllipsoid
    from org.orekit.frames import FramesFactory, TopocentricFrame
    from org.orekit.orbits import CartesianOrbit, OrbitType
    from org.orekit.forces.gravity import J2OnlyPerturbation
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.events import ElevationDetector, EclipseDetector, EventsLogger
    from org.orekit.propagation.events.handlers import ContinueOnEvent
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import IERSConventions, PVCoordinates, Constants
    from org.orekit.utils import TimeStampedPVCoordinates

    utc = TimeScalesFactory.getUTC()
    e = sc.ts._epoch_dt
    epoch = AbsoluteDate(e.year, e.month, e.day, e.hour, e.minute, float(e.second), utc)
    gcrf = FramesFactory.getGCRF()
    cirf = FramesFactory.getCIRF(IERSConventions.IERS_2010, False)
    teoe = cirf.getFrozenFrame(gcrf, epoch, "TEOE")
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
    earth = OneAxisEllipsoid(sc.re * 1e3, sc.f, itrf)
    out = {}

    # --- 0. frame check: our TEOE->GCRS matrix vs Orekit frozen CIRF
    tr = teoe.getStaticTransformTo(gcrf, epoch)
    R_ok = np.array([[tr.transformVector(Vector3D(*e_)).getX(), tr.transformVector(Vector3D(*e_)).getY(),
                      tr.transformVector(Vector3D(*e_)).getZ()] for e_ in np.eye(3)]).T
    from scipy.spatial.transform import Rotation
    out["teoe_frame_diff_rad"] = float(Rotation.from_matrix(R_ok.T @ sc.R_teoe_gcrs).magnitude())

    # --- 1. numerical J2 propagation from the same osculating state
    from ..orbit.cowell import CowellJ2Orbit
    from ..orbit.elements import kep2cart
    d = sc.design
    r0, v0 = kep2cart(d.a_km, d.e, d.i_rad, d.raan_teoe_rad, 0.0, d.u0_rad, sc.mu)
    mine = CowellJ2Orbit(r0, v0, sc.mu, sc.re, sc.j2, sc.t_end + 3600)
    pv0 = PVCoordinates(Vector3D(*(r0 * 1e3)), Vector3D(*(v0 * 1e3)))
    orbit0 = CartesianOrbit(pv0, teoe, epoch, sc.mu * 1e9)
    integ = DormandPrince853Integrator(1e-3, 300.0, 1e-6, 1e-12)
    prop = NumericalPropagator(integ)
    prop.setOrbitType(OrbitType.CARTESIAN)
    prop.setInitialState(SpacecraftState(orbit0))
    prop.addForceModel(J2OnlyPerturbation(sc.mu * 1e9, sc.re * 1e3, sc.j2, teoe))
    # events: elevation detectors for every station + eclipse detector
    loggers = {}
    for st in sc.stations:
        gp = GeodeticPoint(np.radians(st.lat_deg), np.radians(st.lon_deg), st.alt_m)
        topo = TopocentricFrame(earth, gp, st.name)
        det = ElevationDetector(topo).withConstantElevation(np.radians(sc.cfg.constraints.min_elevation_deg)) \
            .withMaxCheck(2.0).withThreshold(1e-4).withHandler(ContinueOnEvent())
        lg = EventsLogger()
        prop.addEventDetector(lg.monitorDetector(det))
        loggers[st.name] = (lg, topo)

    @JImplements("org.orekit.utils.ExtendedPositionProvider")
    class ErfaSun:
        @JOverride
        def getPosition(self, date, frame):
            dt = date.durationFrom(epoch)
            p = sc.sun.position_gcrs_km(np.array([float(dt)]))[0] * 1e3
            return gcrf.getStaticTransformTo(frame, date).transformPosition(Vector3D(*p))

        @JOverride
        def getPVCoordinates(self, date, frame):
            dt = float(date.durationFrom(epoch))
            p = sc.sun.position_gcrs_km(np.array([dt]))[0] * 1e3
            p2 = sc.sun.position_gcrs_km(np.array([dt + 1.0]))[0] * 1e3
            pv = TimeStampedPVCoordinates(date, Vector3D(*p), Vector3D(*(p2 - p)))
            return gcrf.getTransformTo(frame, date).transformPVCoordinates(pv)

    sun_ok = True
    try:
        sun = ErfaSun()
        ecl = EclipseDetector(sun, Constants.SUN_RADIUS, earth).withPenumbra() \
            .withMaxCheck(30.0).withThreshold(1e-4).withHandler(ContinueOnEvent())
        lg_e = EventsLogger()
        prop.addEventDetector(lg_e.monitorDetector(ecl))
    except Exception as exc:  # pragma: no cover
        sun_ok = False
        out["eclipse_detector_error"] = repr(exc)

    t_end = AbsoluteDate(epoch, float(sc.t_end))
    t_start = AbsoluteDate(epoch, float(sc.t_start))
    # sample positions hourly (ephemeris generator) for the propagation comparison
    gen = prop.getEphemerisGenerator()
    prop.propagate(t_end)
    eph = gen.getGeneratedEphemeris()
    ts = np.arange(0.0, sc.t_end, 3600.0)
    dmax = 0.0
    el_err = 0.0
    for tk in ts:
        pv = eph.getPVCoordinates(AbsoluteDate(epoch, float(tk)), teoe)
        p = np.array([pv.getPosition().getX(), pv.getPosition().getY(), pv.getPosition().getZ()]) / 1e3
        rm, _ = mine.state_teoe(np.array([tk]))
        dmax = max(dmax, float(np.linalg.norm(p - rm[0])))
    out["cowell_vs_orekit_max_position_diff_m"] = dmax * 1e3

    # --- 2. topocentric elevation of OUR state, Orekit frames vs ours
    tt = np.linspace(sc.t_start, sc.t_end, 300)
    r_g, _ = sc.state_gcrs(tt)
    for k, st in enumerate(sc.stations):
        topo = loggers[st.name][1]
        e_ok = np.array([np.degrees(topo.getElevation(Vector3D(*(r_g[i] * 1e3)), gcrf, AbsoluteDate(epoch, float(tt[i]))))
                         for i in range(len(tt))])
        e_me = sc.elevation_deg(k, tt)
        el_err = max(el_err, float(np.abs(e_ok - e_me).max()))
    out["elevation_vs_orekit_max_deg"] = el_err

    # --- 3. AOS/LOS events: Orekit vs our event finder on the SAME (Cowell-osculating) trajectory
    from ..scenario import Scenario
    from ..windows import compute_constraints
    sc_c = Scenario(sc.cfg, propagator="cowell_osculating")
    cs_c = compute_constraints(sc_c)
    ev_err = []
    for st in sc.stations:
        lg = loggers[st.name][0]
        times = sorted(float(ev.getState().getDate().durationFrom(epoch)) for ev in lg.getLoggedEvents())
        times = [x for x in times if sc.t_start < x < sc.t_end]
        mine_t = np.sort(sc_c and cs_c.passes[st.name].ravel())
        mine_t = mine_t[(mine_t > sc.t_start + 1e-6) & (mine_t < sc.t_end - 1e-6)]
        if len(times) == len(mine_t):
            ev_err.append(float(np.abs(np.array(times) - mine_t).max()))
        else:
            ev_err.append(float("inf"))
            out[f"event_count_mismatch_{st.name}"] = [len(times), len(mine_t)]
    out["aos_los_vs_orekit_max_s"] = max(ev_err)
    if sun_ok:
        times = sorted(float(ev.getState().getDate().durationFrom(epoch)) for ev in lg_e.getLoggedEvents())
        times = np.array([x for x in times if sc.t_start < x < sc.t_end])
        mine_e = np.sort(cs_c.lit.ravel())
        mine_e = mine_e[(mine_e > sc.t_start + 1e-6) & (mine_e < sc.t_end - 1e-6)]
        if len(times) == len(mine_e):
            out["penumbra_events_vs_orekit_max_s"] = float(np.abs(times - mine_e).max())
        else:
            out["penumbra_event_count"] = [len(times), len(mine_e)]
    ok = (out["cowell_vs_orekit_max_position_diff_m"] < 50.0 and out["elevation_vs_orekit_max_deg"] < 1e-3
          and out["aos_los_vs_orekit_max_s"] < 0.5)
    out["verdict"] = "PASS" if ok else "FAIL"
    return out
