# aule-rpod — Docking-window analysis & RPOD final-approach design

Flight-dynamics pipeline for the Aule Space mission-design assignment: find every epoch in the
10-day mission window at which a 9-minute final approach to a 550 km SSO target satisfies the
illumination, ground-station and eclipse constraints (Part I), then design and justify the
approach from the first window to the 10 m hold point and docking (Part II).

## Answer in one paragraph

**As specified, the mission window contains no valid docking window.** The target's docking
port precesses about +V-bar at exactly 1.2 deg/s with a 10 deg nutation (closed-form law
recovered from the attitude file, residual 2e-11 rad), so the port sun angle rides a ±5–15 deg,
5-minute oscillation on top of the orbital sun ramp. With LTAN 10:30 in July (beta ≈ −25 deg)
the 30–60 deg band is met only on the southern day arc, which only Wellington covers; its best
pass (12 Jul 23:09:58Z, 89 deg max elevation) offers 462 s of the required 540 s. The result is
robust to every modelling choice tested (12 alternative readings, 7,848 orbit-phase/LTAN cases)
and was verified against Orekit (event times to 1 ms) and an independent astropy chain. The
earliest window under the **minimum relaxation** (illumination band widened by 3.41 deg, all
other constraints hard) opens **17 Jul 2026 22:53:08.266 UTC (Wellington, 583.9 s)**. There, a
**+V-bar approach** — the target holding its port on +V-bar — reaches the 10 m hold point and
docks at 23:02:21.9 UTC for **1.184 m/s** (four impulses + 2nL compensation), keeping the Sun at
41–51 deg on the port and ≥ 129 deg off the chaser camera during the critical phase. In 5,000
closed-loop Monte Carlo runs every contact falls inside the translational limits of the IDSS initial-contact conditions;
the recommended delta-v allocation is 1.80 m/s.

## Quick start

```bash
make install        # pip install -e ".[dev,validation]"  (Python >= 3.10)
make test           # 24 unit/regression tests (~10 s)
make selftest       # built-in health check with result codes (~2 s)
make all            # every product + provenance manifest (~16 min, dominated by the 7,848-case maps)
make report         # typeset PDF report (pandoc + XeLaTeX) -> docs/report/report.pdf
```

Individual stages: `aule-rpod part1 | part2 | studies [--quick] | validate`.
Optional V&V data (`make data`): Orekit data directory and JPL DE421 arrays, built from public sources.

## Repository layout

```
config/mission.yaml          single source of truth (every number; ASSUMPTIONs flagged)
data/                        Target_Attitude_Profile.a (+ optional orekit-data/, de421/)
src/aule_rpod/
  config.py                  pydantic schema: typed, range-checked, unknown keys rejected
  timesys.py  frames.py      UTC/TT/UT1 + leap-second guard; IAU 2006/2000A GCRS<->ITRS; STK LVLH
  ephemeris.py               Sun: ERFA epv00 (primary), JPL DE421 + analytic (V&V)
  orbit/                     SSO design, secular-J2 (STK J2Perturbation formulation), Cowell J2,
                             numerically consistent mean-element fit
  attitude/                  STK .a reader, quaternion convention (tested), SLERP/NLERP/law,
                             closed-form precession-law identification
  geometry/                  station look angles, ellipsoidal conical shadow model
  events/                    sample-bracket-refine event finder (1 ms) with hidden-crossing guard;
                             interval algebra
  windows.py relaxation.py   docking-window logic, per-pass diagnostics, minimum relaxations
  sensitivity.py studies.py  feasibility maps, relaxation curves, alternative interpretations
  rpod/                      CW dynamics, waypoint plan builder, strategy trade, safety (missed
                             burns, breakout scan), vectorised closed-loop Monte Carlo (MCC slots,
                             final-translation controller, IDSS contact scoring), nonlinear J2 check
  validation/                cross-checks vs astropy, JPL DE421, Orekit 13
  viz/                       Gantt (static + interactive), analysis figures
  part1.py part2.py cli.py   pipelines and the single entry point
  selftest.py                fast self-test with PASSED/FAILED result codes
tests/                       pytest suite
scripts/                     data preparation (Orekit data, DE421) and submission packaging
docs/report/                 report source (Markdown), PDF build script, typeset report.pdf
outputs/                     generated products (CSV, figures, JSON summaries, manifest)
```

## Architecture

```
mission.yaml ─► MissionConfig (validated) ─► Scenario ─┬─ TimeSystem (UTC/TT/UT1, EOP)
                                                       ├─ EarthOrientation (GCRS↔ITRS)
                                                       ├─ Orbit (secular J2 | Cowell J2)
                                                       ├─ Attitude (SLERP | NLERP | law)
                                                       ├─ Sun (ERFA, Hermite-cached)
                                                       └─ Stations (WGS-84)
          every constraint = continuous g(t) >= 0  ─►  event finder (ms)  ─►  interval algebra
                                                                 │
              Part I: passes ∩ sunlight ∩ sun band ≥ 9 min ──────┤──► CSV / Gantt / summary
              diagnostics, minimum relaxations, studies           │
              Part II: window ─► strategy trade ─► CW plan ─► safety, Monte Carlo, J2 check
```

Design rules: one configuration file; every product timestamped in ISO-8601 UTC and traceable
through `outputs/run_manifest.json` (config fingerprint, source hash, input hash, library
versions); conventions stated once and unit-tested; every model cross-checked against an
independent implementation; seeded randomness; the self-test must pass before products are made.

## Key modelling choices (details and justification in the report)

| Item | Choice |
|---|---|
| Orbit | a = 6928.137 km, e = 0, i = 97.593 deg (J2 sun-synchronous), node at local mean solar time 10:30 at epoch, argument of latitude 0 at epoch (ASSUMPTION) |
| Propagation | first-order secular J2 (STK *J2Perturbation* formulation), frame = true equator & CIO of orbit epoch; Cowell J2 for sensitivity |
| LVLH | STK definition: X radial, Z orbit normal, Y = Z × X (along-track) |
| Quaternions | STK: scalar last, rotation reference → body; port in LVLH = R(q)·p_B |
| Sun / Earth | ERFA epv00 (≈ JPL DE421 to 0.01"); IAU 2006/2000A with measured IERS EOP; WGS-84 |
| Eclipse | conical umbra + penumbra on the WGS-84 ellipsoid; any occultation counts (ASSUMPTION) |
| Stations | altitude 0 m above WGS-84 (ASSUMPTION), 5 deg mask, no refraction |
| Part II | chaser at rest in LVLH at window start (ASSUMPTION), CW dynamics, n from the target state |
| Monte Carlo | navigation 1 % of range (≥ 1 cm) + 1 mm/s, burns 1 % / 0.5 deg; MCC at half of each leg; final translation flown by a 2-s closed-loop controller (20-s lateral time constant, 2 mm/s deadband, 30-s Gauss-Markov navigation errors); contact scored against IDSS IDD Rev. E Table 3.3.1.1-2 |

## Verification & validation (outputs/validation)

| Check | Result |
|---|---|
| Attitude law vs all 14,401 samples | max residual 2.2e-11 rad; SLERP within 0.19 deg of the law |
| Sun: ERFA vs JPL DE421 | 0.01 arcsec |
| Frames vs ERFA direct / astropy | 1e-12 rad / 2 cm, elevation 5e-7 deg |
| Event-finder grid independence (1/2/5 s) | < 1 ms |
| Independent astropy chain, top candidate passes | longest feasible span within 0.7 s |
| Orekit 13 (numerical J2, ElevationDetector, EclipseDetector) | 1.9 cm over 20 days; AOS/LOS 0.9 ms; shadow 14 ms |
| Nonlinear J2 two-spacecraft flight of the CW plan | ≤ 5 cm at every waypoint and at contact |

## Main outputs

| File | Content |
|---|---|
| `outputs/part1/csv/ground_station_passes.csv` | 258 passes: AOS/LOS/TCA, duration, max elevation, azimuths |
| `outputs/part1/csv/eclipse_windows.csv` | 151 eclipses: penumbra/umbra entry & exit |
| `outputs/part1/csv/sun_angle_timeseries_10s.csv` | sun angle on the port every 10 s, band flag, visible Sun fraction |
| `outputs/part1/csv/state_timeline_60s.csv` | full state at the attitude epochs (Sun & port in LVLH, beta, elevations) |
| `outputs/part1/csv/docking_windows.csv` | strict windows (header only: none exist) |
| `outputs/part1/csv/docking_window_candidates.csv` | all 130 qualifying passes ranked, limiting constraints, minimum relaxation |
| `outputs/part1/csv/docking_windows_min_relaxation.csv` | the window used for Part II |
| `outputs/part1/figures/gantt_docking_windows*.{png,pdf,html}` | Gantt chart (static + interactive) |
| `outputs/part2/csv/maneuver_plan.csv` | burns and mid-course-correction slots with UTC epochs and LVLH delta-v components |
| `outputs/part2/csv/strategy_trade.csv` | V-bar / R-bar / H-bar / sun-optimal trade metrics |
| `outputs/part2/analytic_dv_worksheet.json` | every intermediate CW quantity of the analytic delta-v |
| `outputs/sensitivity/` | feasibility maps, relaxation curves, alternative interpretations |
| `outputs/validation/validation_report.{md,json}` | V&V campaign |
| `docs/report/report.pdf` | technical report (assumptions, interpretations, limitations, strategy, delta-v) |
