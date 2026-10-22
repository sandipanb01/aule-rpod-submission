As specified, no epoch in the 10–20 July 2026 mission window supports a 9-minute final approach: the best pass (Wellington, 12 Jul, from 23:09:58 UTC) offers 7.7 of the 9 minutes. The earliest window under the smallest relaxation, an illumination band widened by 3.41°, opens 17 Jul 2026 22:53:08 UTC; there a +V-bar approach docks at 23:02:21.9 UTC for 1.184 m/s.

## 1. Results at a glance

| Item | Result |
|:-------------------------|:-------------------------------------------------------------------------|
| Strict docking windows (all three constraints held for ≥ 9 min) | **None.** 0 of the 130 ground-station passes that last ≥ 9 min contain a 9-min span that is sunlit with the port Sun angle in 30–60° |
| Closest opportunity | Wellington, 12 Jul 2026 23:09:58–23:19:49 UTC, max elevation 89.2°: 462.1 s feasible, 77.9 s short. The start is limited by the Sun angle dropping below 30° |
| Root cause | The port nutates 10° ± 5.0° about +V-bar at 1.2°/s. With β ≈ −25° (LTAN 10:30 in July) the band holds only on the southern day arc, which only Wellington covers |
| Smallest change that creates a window | Widen the illumination band by 3.41° (floor 30° → 26.59°). Changing the elevation mask or the eclipse rule does not help |
| Window used for Part II | Wellington, 17 Jul 2026 22:53:08.266–23:02:52.211 UTC (583.9 s) |
| Approach strategy | +V-bar from ahead; the target holds its port on +V-bar; hold point 10 m ahead of the port |
| Maneuvers (analytic, Clohessy–Wiltshire) | B1 55.69, B2 42.18, B3 12.33, B4 6.00 cm/s + 2.19 cm/s continuous = **1.184 m/s**, plus two mid-course-correction slots |
| Contact | 23:02:21.9 UTC, 30.3 s before loss of signal, at 6 cm/s, with the Sun 41–51° from the port normal and ≥ 129° off the chaser camera |
| Dispersions (5,000 closed-loop runs) | Every run docks inside the IDSS contact envelope (lateral misalignment p99 6.2 cm against 10 cm; closing 5.8–6.2 cm/s against 5–10 cm/s); none enters the keep-out sphere outside the corridor; Δv p99 1.343 m/s |
| Δv allocation | 1.80 m/s: nominal + dispersions + one breakout burn + 20 % margin |
| Verification | 7 independent cross-checks pass, incl. Orekit (pass events within 0.9 ms, shadow events within 14 ms); a nonlinear J2 flight of the plan agrees to ≤ 5 cm |

## 2. Scenario model and assumptions

Given values are used as stated. Every other number is an explicit assumption or is derived, and each assumption is tested for influence in Sections 7 and 9.

| Parameter | Value used | Basis |
|:----------------|:-------------------------------------------|:--------------------------------------|
| Orbit epoch | 1 Jul 2026 06:30:00 UTC; mission window 10 Jul 06:30 to 20 Jul 06:30 UTC | Given |
| Semi-major axis | 6928.137 km (WGS-84 equatorial radius + 550 km) | Given altitude, taken as the mean orbit |
| Eccentricity | 0 | Given |
| Inclination | 97.5930° | Derived: first-order J2 nodal rate equals the mean Sun rate, 0.98564736°/day |
| Right ascension of ascending node | 76.4905° in the true-equator/CIO frame of the orbit epoch | Derived: node at local mean solar time 10:30, i.e. longitude 60.0° E at 06:30 UT1 |
| Argument of latitude at epoch | 0° (target at the ascending node on 1 Jul 06:30 UTC) | **Assumption** (STK orbit-wizard default); swept over 0–360° in Section 7 |
| Propagator | First-order secular J2: Ω, ω and M drift linearly; two-body positions from the drifted elements | "Only J2" as specified; the formulation of STK's J2 Perturbation propagator. Numerical (Cowell) J2 used for sensitivity |
| Earth constants | μ = 398600.4418 km³/s², R = 6378.137 km, f = 1/298.257223563, J2 = 1.08262668×10⁻³ | WGS-84 / EGM96 |
| Time and Earth orientation | UTC with TAI−UTC = 37 s (no leap second in the span); IAU 2006/2000A; measured IERS EOP | Standard |
| Sun | ERFA epv00 geometric geocentric position | Agrees with JPL DE421 to 0.01″ |
| LVLH frame | X radial, Z orbit normal, Y = Z × X (along-track) | STK definition named in the attitude file header |
| Quaternions | [q1, q2, q3, q4], scalar last, rotation from LVLH to body; the port in LVLH is R(q) applied to the port in body axes | STK attitude-file format |
| Attitude between 60-s samples | Spherical linear interpolation (SLERP) | Checked against the exact identified law and STK's linear interpolation: ≤ 0.22° difference |
| Docking-port normal | [0.9962, 0, 0.0872] normalised (input norm 1.0000091); 5.0025° from body +X toward body +Z | Given |
| Sun angle | Angle between the port normal and the target-to-Sun vector, 0–180° | Given |
| Eclipse | Any umbra or penumbra, conical shadow on the WGS-84 ellipsoid | **Assumption** (conservative); shown not to be binding |
| Ground stations | Coordinates as given; height 0 m; 5° elevation mask; no refraction | **Assumption** for height and refraction |
| Docking window | Start epochs t₀ such that [t₀, t₀ + 540 s] lies inside one pass of ≥ 540 s, in full sunlight, with the Sun angle in [30°, 60°] throughout | Interpretation of the assignment text |
| Part II initial state | ρ = (−80, +120, +20) m in LVLH (below, ahead, +orbit-normal side); ρ̇ = 0 | Distances given; **assumption** for the cross-track side (either side gives the same Δv) and for zero relative velocity (natural-drift case in Section 12) |
| Part II target attitude | Port held along the chosen approach axis, fixed in LVLH | Given: "the target can change docking port orientation" |

## 3. What the attitude file encodes

All 14,401 quaternions follow one closed-form law to a maximum residual of 2.2×10⁻¹¹ rad. The target is not tumbling freely: its docking port nutates about +V-bar. The law was fitted by Levenberg–Marquardt on geodesic residuals, seeded from the data's own per-step rotation axes. The attitude is a spin about body +X followed by a precession about LVLH +Y:

```latex
R_{L\leftarrow B}(t) = \mathrm{Rot}\big(\hat{y}_L,\ \omega_p\,(t-t_0)\big)\; R_0\; \mathrm{Rot}\big(\hat{x}_B,\ \omega_s\,(t-t_0)\big)
```

| Parameter | Identified value |
|:------------------------------|:-----------------------------------------------|
| Precession axis, fixed in LVLH | +Y (along-track, +V-bar), to 10⁻¹⁵ |
| Precession rate | 1.200000000°/s (period 300.000 s) |
| Spin axis, fixed in the body | body +X, to 10⁻¹⁵ |
| Spin rate | 0.109384178°/s (period 3291.15 s) |
| Nutation: angle between LVLH +Y and the spin axis | 10.0000000° |

- The port normal sits 5.0025° from body X, so it always lies 10° ± 5.0025° from +V-bar. Its azimuth about +V-bar sweeps once every 5 minutes.
- The Sun angle on the port is therefore the Sun-to-V-bar angle plus a 5-minute oscillation of up to ±15°. Over 9 minutes that oscillation alone spans more than half of the 30°-wide band.
- Between the 60-s samples, SLERP deviates from the exact law by ≤ 0.19° in port direction, and STK's linear interpolation by ≤ 0.22°. Window lengths change by < 1 s.
- Convention check: reading the quaternions as body-to-reference makes the port sweep the whole Y–Z plane every 5 minutes, and the 30–60° band then never lasts more than 86 s. Together with the STK file specification, this confirms the reference-to-body reading.

![Docking-port normal seen along +V-bar over the first hour (left) and its offset from +V-bar over three hours (right).](FIG:attitude)

## 4. Method

Every constraint is a continuous function g(t) that is non-negative exactly when the constraint holds. Events are located to 1 ms and the satisfied sets are combined with interval algebra.

```latex
\begin{aligned}
g_{\mathrm{vis},k}(t) &= \mathrm{el}_k(t) - 5^\circ \\
g_{\mathrm{lit}}(t) &= d(t) - \big(\alpha_\oplus(t) + \alpha_\odot(t)\big) \\
g_{30}(t) &= \theta(t) - 30^\circ, \qquad g_{60}(t) = 60^\circ - \theta(t), \qquad \theta(t) = \arccos\big(\hat{p}_L(t)\cdot\hat{s}_L(t)\big)
\end{aligned}
```

- The elevation el is measured above station k's WGS-84 local horizontal. The target position is rotated from GCRS to ITRS with the IAU 2006/2000A CIO-based chain and measured Earth-orientation parameters.
- Shadow, conical model: α⊕ = asin(R⊕/|r|) and α⊙ = asin(R⊙/|s − r|) are the apparent radii of Earth and Sun seen from the target, and d is the angle between their centres. The lighting function is non-negative only when no part of the solar disc is hidden. Oblateness is handled by scaling z by a/b, which maps the WGS-84 ellipsoid onto a sphere and keeps rays straight.
- p̂ is the port normal rotated into LVLH by the attitude R(t), and ŝ the target-to-Sun unit vector in LVLH.

Each g is sampled every 2 s. Every sign change is bracketed and refined by vectorised bisection to 1 ms, and a parabolic test on each triple of same-sign samples catches double crossings hidden inside one step. Results agree to < 1 ms on 1-s, 2-s and 5-s grids.

With P(k, j) the j-th pass of station k, L the full-sunlight set, B the Sun-angle band and D = 540 s, the docking windows are:

```latex
F_{k,j} = P_{k,j}\cap L\cap B,\qquad
\mathcal{W}=\Big\{\,[a,\;b-D]\ :\ [a,b]\in F_{k,j},\ \ b-a\ge D,\ \ |P_{k,j}|\ge D\,\Big\}
```

Each element of W is an interval of admissible start epochs. For every pass the pipeline also reports the longest feasible span, its deficit to D, the constraint that bounds each end, and the smallest relaxation that would close the gap (sliding-window extremes of θ on a 0.5-s grid).

## 5. Part I results

| Constraint | 10-day statistics |
|:------------------------------|:------------------------------------------------------------|
| Svalbard passes | 140; 80 last ≥ 9 min; longest 9.98 min |
| Bengaluru passes | 32; 15 last ≥ 9 min; longest 9.75 min |
| Wellington passes | 45; 19 last ≥ 9 min; longest 9.84 min |
| Los Angeles passes | 41; 16 last ≥ 9 min; longest 9.80 min |
| Eclipses | 151 (the last cut short by the end of the window); every full eclipse lasts 34.2–34.3 min, of which 33.9 min umbra; β from −24.9° to −25.2° |
| Sun angle within 30–60° | 23.2 % of the time; 1,001 intervals, 63 of them ≥ 9 min |
| Full sunlight and band together for ≥ 9 min | 29 intervals |
| All constraints for ≥ 9 min within one pass | **0** |

A pass lasts ≥ 9 min only if its maximum elevation exceeds about 28°. The longest feasible span inside any pass is 462 s at Wellington and 49 s at Bengaluru. At Svalbard and Los Angeles it is 0 s: during their passes the port Sun angle never enters the band.

![Constraint Gantt chart over the 10-day mission window. Passes shorter than 9 min are faint; orange marks the longest feasible span of the six closest passes; the diamond is the window created by the minimum relaxation. An interactive, zoomable version is delivered with the outputs.](FIG:gantt)

| # | Pass AOS (UTC, Wellington) | Max el. | Longest feasible | Short by | Bounded by (start / end) | Band widening needed |
|:--|:-------------|:------|:--------|:--------|:----------------|:------|
| 1 | 12 Jul 23:09:58 | 89.2° | 462.1 s | 77.9 s | θ < 30° / LOS | 7.64° |
| 2 | 11 Jul 23:13:21 | 81.5° | 417.0 s | 123.0 s | θ < 30° / LOS | 8.06° |
| 3 | 17 Jul 22:53:08 | 55.8° | 387.3 s | 152.7 s | θ < 30° / θ > 60° | **3.40°** |
| 4 | 16 Jul 22:56:29 | 61.8° | 337.7 s | 202.3 s | θ < 30° / LOS | 3.69° |
| 5 | 10 Jul 23:16:44 | 74.0° | 333.3 s | 206.7 s | θ < 30° / LOS | 9.59° |
| 6 | 13 Jul 23:06:36 | 83.1° | 298.0 s | 242.0 s | θ < 30° / LOS | 7.55° |

All 130 qualifying passes, with their limiting constraints, are ranked in docking_window_candidates.csv.

![Anatomy of the six closest passes: the port Sun angle (black) is the cone-axis angle (grey ramp) plus the 5-minute nutation; in every case it dips below 30° early in the pass.](FIG:anatomy)

## 6. Why no window exists

- In LVLH the Sun stays about 25° out of the orbit plane (β ≈ −25°) while its in-plane component turns once per orbit. The angle between the Sun and +V-bar, the cone axis of the port, lies inside 30–60° only on the ascending (day) arc between about 51° S and 12° S.
- A 9-minute pass needs the target within about 7.5° of a station. Only Wellington, at 41.2° S, puts such passes on that arc, between 49° S and 34° S. There the cone-axis angle is only 32–42°, so the ±5–15° nutation drives θ below 30° in the first minutes of every pass.
- Svalbard, Bengaluru and Los Angeles see the day arc where the Sun is 76–148° from +V-bar. No orbit phase changes this: the Sun–V-bar geometry is tied to local time, and the nutation phase to absolute time.

![Sun angles on the ascending day arc versus latitude, with each station's coverage below. The band is reachable only where Wellington alone has ≥ 9-min passes.](FIG:geometry)

| Relaxation (all else at specification) | Windows in 10 days | First window (UTC) |
|:------------------------------------------|:--------|:------------------------------|
| Band widened by 3.41° (floor 26.59°, ceiling 63.41°) | 1 | 17 Jul 22:53:08, Wellington |
| Floor lowered to 26° | 2 | 16 Jul 22:56:29, Wellington |
| Floor lowered to 22° | 4 | 12 Jul 23:09:58, Wellington |
| Required approach shortened to 7.5 min | 1 | 12 Jul 23:12:07, Wellington |
| Required approach shortened to 6.0 min | 3 | 11 Jul 23:16:14, Wellington |
| Elevation mask lowered to 0° | 0 | none |
| Penumbra allowed during the approach | 0 | none |

The smallest uniform widening that opens a window is 3.403°, at which the whole 17 July pass becomes feasible; Part II uses 3.41°.

Recommendation: relax the illumination band, not the safety constraints. The ground-station link carries the go/no-go and abort commands, full sunlight protects power and sensors, and a shorter approach means higher closing speeds. The band edges are about lighting quality: at the relaxed edges (26.6° and 63.4°) the chaser's shadow at 10 m still falls 5 m from the port, and the face still receives 45 % of normal-incidence light. Part II therefore uses the window created by the minimum band widening; with the target re-oriented, the port then sees the Sun at 41–51°, inside the original band.

![Number of strict windows as each requirement alone is relaxed.](FIG:relaxation)

## 7. Robustness to the unspecified orbit parameters

The no-window result holds for 71 of 72 assumed orbit phases and for every tested reading of the orbit and attitude specification.

- Initial argument of latitude: at LTAN 10:30 only 1 of 72 phases (5° steps) yields a window. The longest feasible span ranges from 391 s to 549 s across phases.
- LTAN: with LTAN between 09:10 and 10:00, 81–99 % of orbit phases give windows; the share falls to 26 % at 10:18 and to 1 of 72 at 10:30, and is zero from 10:36. Over a full day, most phases give windows only for LTAN of about 01:45–03:30, 09:00–10:10, 15:00–16:30 and 19:45–22:30.

![Strict-window count versus LTAN and initial argument of latitude around the specified orbit (circle).](FIG:mapzoom)

| Reading of the specification | Strict windows | Longest feasible span |
|:----------------------------------------------------|:--------------|:----------|
| As specified (baseline) | 0 | 462.1 s |
| Descending node at 10:30 instead of ascending node | 3 (Svalbard; first 11 Jul 09:30:07) | 560.4 s |
| Quaternions read body → reference | 0 | 86.0 s |
| Attitude interpolated linearly (STK default) | 0 | 461.9 s |
| Attitude from the exact identified law | 0 | 461.9 s |
| Numerical J2, entered elements taken as osculating | 0 | 455.6 s |
| Numerical J2, entered elements taken as mean | 0 | 460.6 s |
| Penumbra allowed | 0 | 462.1 s |
| LTAN referenced to the apparent Sun | 0 | 461.2 s |
| Argument of latitude at epoch 90° / 180° / 270° | 0 / 0 / 0 | 442.8 / 442.9 / 504.8 s |

If a reference solution reports windows at Svalbard near 09:30 UTC, it corresponds to a descending node at 10:30; those windows are listed in alternative_interpretations_windows.csv.

## 8. Verification and validation

Every model was checked against an independent implementation. Tolerances come from what the analysis needs (event times to about 1 s, angles to about 0.01°).

| Check | Reference | Result |
|:----------------------|:----------------------------|:----------------------------------|
| Attitude law and quaternion convention | All 14,401 samples; canonical rotation test | Residual 2.2×10⁻¹¹ rad; SLERP within 0.19° of the law |
| Sun ephemeris | JPL DE421; analytic almanac model | 0.01″ (DE421); 19″ (almanac) |
| Earth orientation | ERFA direct; astropy | 1×10⁻¹² rad; 2 cm in position; 5×10⁻⁷° in elevation |
| Orbit model | Numerical J2 (mean and osculating readings); LTAN at all 285 nodes | Flown LTAN within 0.15 s of 10:30:00 at every node; mean-reading phase drift 0.024°/day |
| Event finder | 1-s, 2-s and 5-s grids | Identical to < 1 ms |
| End-to-end, top candidate passes | Independent astropy chain (frames, Sun) with the exact attitude law | Longest feasible span within 0.71 s |
| Propagation, visibility and shadow | Orekit 13 (numerical J2, ElevationDetector, EclipseDetector) | 1.9 cm over 20 days; AOS/LOS within 0.9 ms; shadow within 14 ms |
| Part II dynamics | Two-spacecraft numerical J2 flight of the plan | ≤ 5.1 cm at every waypoint and at contact |

The code also carries 24 unit and regression tests and a built-in self-test with result codes that must pass before any product is written.

## 9. Limitations audit

| Limitation | Effect on the result | Status |
|:---------------------------|:----------------------------------|:------------------------------------------|
| Orbit phase at epoch not given | Moves every pass and eclipse in time | Bounded: only 1 of 72 phases (50°) gives a window, so the conclusion holds for 71 of 72 (Section 7) |
| "550 km circular" read as mean or osculating | The osculating reading puts the target 99° further along-track at mission start | Bounded: 0 windows under both readings (Section 7) |
| Only J2, as specified (no drag, radiation pressure or third bodies) | Drag at 550 km would shift along-track timing over 19 days | Bounded by the phase sweep; in operations a fresh orbit determination before the window replaces a long propagation [Vallado 2013] |
| Attitude known only every 60 s | Interpolation choice matters between samples | Quantified: ≤ 0.19° (SLERP) and ≤ 0.22° (linear) from the exact law; spans change < 1 s |
| Earth-shadow model | A spherical Earth moves shadow edges by 6.8 s | Resolved with the WGS-84 ellipsoid (Orekit agreement 14 ms); atmospheric refraction not modelled [Montenbruck & Gill 2000]; eclipse never binds |
| Station height 0 m, no refraction | Refraction near 5° elevation lengthens passes by a few seconds | Not modelled; lowering the mask to 0° still gives no window |
| Geometric Sun direction | Aberration ≤ 20.5″ | Negligible |
| Constraint interpretation (penumbra counted as eclipse) | Could move window edges | Tested: allowing penumbra still gives 0 windows |
| Part II: linear Clohessy–Wiltshire dynamics | Errors grow with range²/radius | Quantified: ≤ 5.1 cm against a nonlinear J2 flight |
| Part II: point-mass target with the port at its centre | A real port offset shifts the hold point by a constant | Negligible Δv effect |
| Part II: approach squeezed into a 9.7-min window | Closing speed 0.52 m/s near 100 m range, faster than crewed-vehicle practice [Fehse 2003] | Mitigated: passively safe far leg, slow corridor (0.12 m/s), breakout burn verified from every point |
| Part II: stop burn B3 at 10 m | Exhaust points 7.6° from the target (plume impingement) | Mitigation: fire B3 on a canted thruster pair (≥ 25°) or brake before 10 m; needs the thruster layout |
| Part II: navigation, burn-error and controller models | Drive the dispersions and the closed-loop Δv | Assumed: navigation 1 % of range (≥ 1 cm) and 1 mm/s, burns 1 % and 0.5°, controller gains and deadband. All 5,000 runs dock inside the IDSS translational envelope |
| Part II: chaser attitude not modelled | The IDSS angular conditions (4°, 0.2°/s) are not assessed | Out of scope here; needs the chaser's attitude-control design |
| Part II: initial relative velocity not given | Changes B1 | Quantified: a natural-drift start adds 0.082 m/s (Section 12) |

## 10. Part II: approach strategy

+V-bar is selected. Only two candidate axes keep the port lit inside the band, the Sun behind the camera, the Earth out of view and the direct path outside the keep-out sphere: +V-bar and an inclined Sun-optimal axis. Of the two, +V-bar costs 19 % less Δv and its hold point is a natural equilibrium.

During the selected window the Sun direction in LVLH moves from (0.173, 0.889, −0.425) to (0.669, 0.610, −0.425). The Sun is ahead of the target, rises from 10° to 42° above the local horizontal, and sits 25° to the −H side.

All candidates fly the same timeline: two Clohessy–Wiltshire legs (235 s and 122 s) via a corridor-entry point 25 m out, a 30-s hold at 10 m, and a 10-m final translation at 6 cm/s.

| Approach axis (port normal) | Δv [m/s] | Sun on port, hold + final | Sun off camera | Camera background | Miss if B2 missed | Path enters keep-out sphere | Verdict |
|:-------------------|:------|:------------|:--------|:----------|:--------|:---------|:--------------------|
| +V-bar, from ahead | 1.184 | 40.9–50.8° | ≥ 129° | deep space | 10.0 m | no | **selected** |
| Sun-optimal axis (0.12, 0.36, −0.93) | 1.458 | 43.9–46.1° | ≥ 134° | deep space | 25.0 m | no | feasible; +23 % Δv; off-axis station keeping |
| +R-bar, from below | 1.374 | 120–130° (dark) | ≥ 50° | deep space | 24.8 m | no | rejected: port in shadow |
| H-bar, from +H | 1.453 | 115° (dark) | ≥ 65° | deep space | 25.0 m | no | rejected: port in shadow |
| H-bar, from −H | 1.539 | 64.9° | ≥ 115° | deep space | 25.0 m | yes | rejected: out of band |
| −R-bar, from above | 1.649 | 49.6–60.1° | ≥ 120° | Earth | 25.0 m | yes | rejected: Earth background |
| −V-bar, from behind | 1.772 | 129–139° (dark) | ≥ 41° | deep space | 25.0 m | yes | rejected: port in shadow |

- **Trajectory safety.** V-bar points are natural equilibria: the 10-m hold costs no station keeping and does not drift. If the corridor-entry burn B2 is missed, the chaser still misses the target by 10 m. A single 15 cm/s posigrade breakout burn is safe from every point of the corridor, hold and final translation (Section 12).
- **Δv.** The chaser already sits ahead of the target, so +V-bar gives the shortest transfer: 1.184 m/s, while the alternatives cost 16–50 % more.
- **Illumination.** With the port on +V-bar, θ climbs from 27° to 52° over the window and stays at 41–51° within 10 m. The Sun is ≥ 129° off the camera boresight, so there is no glare, and at θ ≥ 41° the chaser's shadow falls ≥ 8.7 m from the port. The camera looks along the horizon, above the Earth's limb, at deep space.
- **Precedent.** The architecture matches crewed ISS visiting vehicles. Crew Dragon's first flight approached the forward port along +V-bar through waypoints at 150 m and 20 m, with a back-away test, a go/no-go poll and a hold for better lighting ([Planetary Society](https://www.planetary.org/articles/crew-dragon-docks)). MEV-1 used 80 m, 20 m and 1 m holds with ground go/no-go ([NASASpaceflight](https://www.nasaspaceflight.com/2020/02/northrop-grumman-history-mission-extension-vehicle-docks-satellite/)).

![Sun angle on the re-oriented port for each candidate axis over the window, selected axis first; green marks the hold and final translation.](FIG:trade)

## 11. Part II: maneuver sequence and analytic Δv

The plan flies P0 (−80, 120, 20) m → WP1 (0, 25, 0) m on the 25-m keep-out sphere → HP (0, 10, 0) m → contact at the port. Relative motion follows the Clohessy–Wiltshire equations, with x radial, y along-track and z cross-track:

```latex
\ddot{x} = 3n^2x + 2n\dot{y} + u_x,\qquad \ddot{y} = -2n\dot{x} + u_y,\qquad \ddot{z} = -n^2 z + u_z
```

With c = cos nt and s = sin nt, the state-transition matrix blocks are:

```latex
\begin{aligned}
\Phi_{rr}&=\begin{bmatrix}4-3c&0&0\\ 6(s-nt)&1&0\\ 0&0&c\end{bmatrix}, &
\Phi_{rv}&=\begin{bmatrix}s/n&2(1-c)/n&0\\ -2(1-c)/n&(4s-3nt)/n&0\\ 0&0&s/n\end{bmatrix},\\[4pt]
\Phi_{vr}&=\begin{bmatrix}3ns&0&0\\ -6n(1-c)&0&0\\ 0&0&-ns\end{bmatrix}, &
\Phi_{vv}&=\begin{bmatrix}c&2s&0\\ -2s&4c-3&0\\ 0&0&c\end{bmatrix}
\end{aligned}
```

Each leg is a fixed-time two-impulse transfer:

```latex
\begin{aligned}
\mathbf{v}_0^{+}&=\Phi_{rv}^{-1}(T)\big[\mathbf{r}_f-\Phi_{rr}(T)\,\mathbf{r}_0\big], &
\Delta\mathbf{v}_1&=\mathbf{v}_0^{+}-\mathbf{v}_0^{-},\\
\mathbf{v}_f^{-}&=\Phi_{vr}(T)\,\mathbf{r}_0+\Phi_{vv}(T)\,\mathbf{v}_0^{+}, &
\Delta\mathbf{v}_2&=\mathbf{v}_f^{+}-\mathbf{v}_f^{-}
\end{aligned}
```

Holding a straight V-bar translation at speed v requires a radial acceleration of 2nv against the Coriolis term, so over a length L:

```latex
\Delta v_{\mathrm{cont}}=\int_0^{L/v}2nv\,\mathrm{d}t = 2nL
```

Worked numbers. The mean motion of the target's LVLH frame at the window start is n = |r × v|/|r|² = 1.0948237×10⁻³ rad/s (period 95.65 min).

- **Leg 1** (T₁ = 235 s, nT₁ = 0.257284 rad): the free-drift image of the start point is Φrr r₀ = (−87.900, 121.358, 19.342) m, so the transfer must supply (87.900, −96.358, −19.342) m. Solving with Φrv gives v₀⁺ = (0.457484, −0.306458, −0.083220) m/s, which is Δv₁ because the chaser starts at rest. It arrives at WP1 with v⁻ = (0.219607, −0.481630, −0.086053) m/s.
- **Leg 2** (T₂ = 122 s, nT₂ = 0.133568 rad): from WP1 to HP the departure velocity is (0.016349, −0.122223, 0) m/s, so Δv₂ = (−0.203258, 0.359407, 0.086053) m/s. Arrival at HP is at (−0.016349, −0.122223, 0) m/s, and Δv₃ cancels it.
- **Final translation**: Δv₄ = 0.06 m/s toward the port, then 2nv = 1.314×10⁻⁴ m/s² radially for 166.7 s, which is 2nL = 0.021896 m/s for L = 10 m.

| Burn | Epoch (UTC, 17 Jul) | t [s] | Δv radial | Δv along-track | Δv cross-track | Δv magnitude [m/s] | Purpose |
|:------|:--------------|:-------|:----------|:----------|:----------|:-------------|:------------------------|
| B1 | 22:53:08.266 | 0 | +0.457484 | −0.306458 | −0.083220 | 0.556897 | leave P0; null the 20-m cross-track offset |
| MCC1 | 22:55:05.766 | 117.5 | 0 | 0 | 0 | 0 (p99 0.059) | re-target WP1 from the navigated state |
| B2 | 22:57:03.266 | 235 | −0.203258 | +0.359407 | +0.086053 | 0.421773 | brake at WP1 and start the 12 cm/s corridor leg |
| MCC2 | 22:58:04.266 | 296 | 0 | 0 | 0 | 0 (p99 0.030) | re-target the hold point |
| B3 | 22:59:05.266 | 357 | +0.016349 | +0.122223 | 0 | 0.123311 | stop at the 10-m hold point |
| hold | 22:59:05–22:59:35 | 357–387 | 0 | 0 | 0 | 0 | go/no-go poll; V-bar equilibrium |
| B4 | 22:59:35.266 | 387 | 0 | −0.060000 | 0 | 0.060000 | start the final translation at 6 cm/s |
| C | 22:59:35–23:02:21.9 | 387–553.7 | +1.314×10⁻⁴ m/s² | 0 | 0 | 0.021896 | radial compensation (2nL) |
| contact | 23:02:21.932 | 553.7 | | | | | 6 cm/s, 30.3 s before loss of signal |

| Δv budget item | Δv [m/s] |
|:-------------------------------------------------------|:------|
| Nominal plan (B1–B4 + compensation) | 1.184 |
| Dispersions: mid-course corrections and closed-loop final translation (Monte Carlo p99 − nominal) | 0.159 |
| Breakout reserve (one 15 cm/s burn) | 0.150 |
| Subtotal | 1.493 |
| System margin, 20 % | 0.299 |
| **Recommended allocation** | **1.80** |

The split between the two legs was chosen by sweeping T₁ with T₁ + T₂ fixed by the window: longer first legs save Δv but enlarge the stop burn B3 and shrink the miss distance if B2 is missed. The chosen point minimises Δv subject to B3 ≤ 15 cm/s and a missed-B2 miss distance ≥ 10 m.

![Approach in the orbital plane with the near field enlarged, the cross-track view and the range history. Grey curves show the free drift if B2 or B3 were missed.](FIG:trajectory)

![Timeline inside the docking window and the Sun angle on the port.](FIG:timeline)

![Choice of the transfer-time split: total Δv, the stop burn B3 and the closest approach if B2 is missed. The miss distance is the binding limit.](FIG:tltrade)

## 12. Part II: safety, dispersions and robustness

| Failure (single missed burn) | Range at failure | Closest approach in free drift (2 orbits) | Outcome |
|:-------------------------|:-----------|:-----------------------------------|:---------------------|
| B1 not executed | 145.6 m | 145.6 m; drifts away below and ahead | safe |
| B2 not executed (corridor entry) | 25.0 m | 10.0 m after 44 s | passively safe |
| B3 not executed (hold stop) | 10.0 m | 2.2 m after 78 s | breakout burn required |
| B4 not executed | 10.0 m | stays at 10 m (V-bar equilibrium) | safe |

A breakout burn of 0.15 m/s was simulated from every 2-s point of the corridor, hold and final translation (down to 0.3 m range) in each candidate direction, and propagated for two orbits:

| Breakout direction | Worst further closing after the burn | Closest return, 10 min to 2 orbits after | Safe |
|:---------------------------|:-----|:-----|:----|
| +Y posigrade (back along +V-bar) | 0.0 m | 28.3 m | yes |
| Along the line of sight, away from the target | 0.0 m | 27.8 m | yes |
| ±X radial | 5.1–6.3 m | 10.0 m | no |
| ±Z out of plane | 5.7 m | 10.0 m | no |
| −Y retrograde | 24.0 m | 79.3 m | no |

The +Y posigrade burn is the designated breakout: it never closes the range and sends the chaser above and behind the target. Radial and out-of-plane burns give closed relative orbits that return to the trigger point after one orbit, and a retrograde burn first closes a further 24 m.

Each Monte Carlo run re-targets every burn from a navigated state (1 % of range per axis, at least 1 cm, and 1 mm/s) and executes it with 1 % magnitude and 0.5° pointing errors (1σ); both legs carry a mid-course correction at half time. The 10-m final translation is flown closed-loop: every 2 s the controller commands 6 cm/s along the corridor plus a lateral correction with a 20-s time constant, and fires only when the change exceeds 2 mm/s, while the navigation errors wander with a 30-s correlation time. Contact conditions are scored against the translational limits of the IDSS initial-contact table.

| Monte Carlo quantity (5,000 closed-loop runs) | Mean | p99 | Limit |
|:-----------------------------------|:-----------|:---------------|:-------------------------|
| Arrival error at the 10-m hold point | 0.30 m | 0.64 m | |
| Residual velocity after the stop burn B3 | | 4.8 mm/s | |
| Mid-course corrections, leg 1 / leg 2 | 27 / 13 mm/s | 59 / 30 mm/s | |
| Closest approach on leg 1 | | 23.0 m (p1) | enters the 25-m sphere only inside the corridor |
| Runs entering the keep-out sphere outside the ±10° corridor | 0 % | | 0 % |
| Inside the corridor at the end of the hold and throughout the final translation | 100 % | | |
| Final-translation Δv (median 30 corrections) | 0.145 m/s | 0.170 m/s | nominal 0.082 m/s |
| Closing speed at contact | 6.0 cm/s | 5.77–6.24 cm/s (p1–p99) | IDSS 5–10 cm/s |
| Lateral rate at contact | | 2.8 mm/s | IDSS 40 mm/s |
| Lateral misalignment at contact | 2.6 cm | 6.2 cm | IDSS 10 cm |
| Runs inside the IDSS translational contact envelope | 100 % | | |
| Total Δv | 1.289 m/s | 1.343 m/s | nominal 1.184 m/s |

The final translation takes 158–175 s (p1–p99) instead of the nominal 166.7 s, so contact stays at least 22 s before loss of signal.

![Monte Carlo distributions: arrival error at the hold point, lateral misalignment and closing speed at contact against the IDSS limits, and total Δv.](FIG:mc)

The same plan flown with both spacecraft integrated numerically under J2 differs from the CW prediction by 1.6 cm at WP1, 2.9 cm at the hold point and 5.1 cm at contact (0.12 mm/s in velocity).

| Window | Available | Best axis | Δv (+V-bar) | Sun on port, hold + final |
|:-----------------------------------------|:----------|:---------|:----------|:--------------|
| 17 Jul 22:53:08, minimum relaxation (used) | 583.9 s | +V-bar | 1.184 m/s | 40.9–50.8° |
| 16 Jul 22:56:29, floor 26° | 586.6 s | +V-bar | 1.188 m/s | 40.7–50.6° |
| 12 Jul 23:12:07, longest strictly feasible span | 462.1 s | +V-bar | 1.990 m/s | 39.8–49.6° |

The choice of +V-bar does not depend on which candidate window is used. If the chaser instead started on a natural circular orbit 80 m below the target, drifting ahead at 0.131 m/s, B1 would rise to 0.639 m/s and the total to 1.266 m/s.

## 13. Deliverables

| File | Content |
|:------------------------------------------------|:------------------------------------------------|
| src/aule_rpod/, config/mission.yaml, tests/, Makefile, README.md | Python source, single configuration file, 24 tests, one-command build |
| outputs/part1/csv/ground_station_passes.csv | 258 passes: AOS, LOS, time of closest approach, duration, maximum elevation, azimuths |
| outputs/part1/csv/eclipse_windows.csv | 151 eclipses: penumbra and umbra entry and exit times |
| outputs/part1/csv/sun_angle_timeseries_10s.csv | Port Sun angle every 10 s with band and sunlight flags |
| outputs/part1/csv/state_timeline_60s.csv | Full state at the attitude epochs: Sun and port in LVLH, β, elevations |
| outputs/part1/csv/docking_windows.csv | Strict windows (header only: none exist) |
| outputs/part1/csv/docking_window_candidates.csv | All 130 qualifying passes ranked, with limiting constraints and minimum relaxations |
| outputs/part1/csv/docking_windows_min_relaxation.csv | The window used in Part II |
| outputs/part1/figures/gantt_docking_windows (PNG, PDF, interactive HTML) | Gantt chart |
| outputs/part2/csv/maneuver_plan.csv | Burns and correction slots with UTC epochs and LVLH Δv components |
| outputs/part2/csv/strategy_trade.csv | Trade metrics for every approach axis |
| outputs/part2/csv/ (other files) | 1-s trajectory with Sun angles, missed-burn and breakout scans, plume geometry, nonlinear J2 flight, window robustness |
| outputs/part2/analytic_dv_worksheet.json | Every intermediate quantity of the analytic Δv |
| outputs/sensitivity/, outputs/validation/ | Feasibility maps, relaxations, alternative readings, V&V reports |
| outputs/run_manifest.json | Provenance of the run: configuration fingerprint, source and input hashes, library versions |
| docs/report/ | This report as PDF and its Markdown source |

## Sources

- Ansys STK Help: [vehicle local coordinate axes](https://help.agi.com/stk/Content/gator/eq-coordsys.htm) and [trajectory axes types](https://help.agi.com/stkdevkit/Content/DocX/AgSTKVgtLib~Enumerations~AgECrdnTrajectoryAxesType_EN.html) (LVLH definition); [attitude file format](https://help.agi.com/stk/Content/stk/importfiles-01.htm) (quaternion convention); [J2/J4 propagators](https://help.agi.com/stk/Content/stk/vehSat_orbitProp_2bodyJ2J4.htm) and [ellipse options](https://help.agi.com/stk/Content/stk/vehSat_orbitProp_J2J4EllipMotion.htm).
- [Crew Dragon Demo-1 docking](https://www.planetary.org/articles/crew-dragon-docks), The Planetary Society, 2019.
- [MEV-1 docks with Intelsat 901](https://www.nasaspaceflight.com/2020/02/northrop-grumman-history-mission-extension-vehicle-docks-satellite/), NASASpaceflight, 2020.
- [International Docking System Standard (IDSS) Interface Definition Document, Revision E](https://www.internationaldockingstandard.com/download/IDSS_IDD_Revision_E_TAGGED.pdf), 2016, Table 3.3.1.1-2 (initial contact conditions).
- Clohessy, W. H. and Wiltshire, R. S., Terminal guidance system for satellite rendezvous, Journal of the Aerospace Sciences 27(9), 1960.
- Fehse, W., Automated Rendezvous and Docking of Spacecraft, Cambridge University Press, 2003.
- Vallado, D. A., Fundamentals of Astrodynamics and Applications, 4th ed., Microcosm Press, 2013.
- Montenbruck, O. and Gill, E., Satellite Orbits, Springer, 2000.
- Markley, F. L. and Crassidis, J. L., Fundamentals of Spacecraft Attitude Determination and Control, Springer, 2014.
- Folkner, W. M. et al., The Planetary and Lunar Ephemeris DE421, IPN Progress Report 42-178, 2009.
- Orekit 13 flight-dynamics library; IAU SOFA/ERFA; IERS Earth-orientation data.
