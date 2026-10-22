# Aule RPOD — Docking-Window Analysis and Terminal-Approach Design

A reproducible astrodynamics analysis of rendezvous, proximity operations, and docking (RPOD). The project evaluates whether a specified orbital scenario provides a continuous docking opportunity under geometric, illumination, eclipse, and ground-station constraints. When the baseline constraints cannot be satisfied, it quantifies the feasibility gap and evaluates a conditional terminal-approach design.

## Key result

**The baseline scenario produces no valid 540-second docking window in the 10-day search.**

The longest strictly feasible candidate span is 462.1 seconds, leaving a 77.9-second deficit. Under a tested relaxation that widens the solar-port-angle band by 3.41°, the analysis identifies a 583.9-second candidate window at Wellington on 17 July 2026.

This distinction matters: the relaxed window is a conditional engineering result, not proof that the original mission constraints are feasible.

| Analysis result | Value |
|---|---:|
| Baseline search horizon | 10 days |
| Required continuous docking window | 540 s |
| Valid baseline windows | 0 |
| Longest strict candidate span | 462.1 s |
| Shortfall against requirement | 77.9 s |
| Minimum tested illumination-band relaxation | 3.41° |
| Candidate span after relaxation | 583.9 s |
| Designed approach direction | +V-bar |
| Nominal modeled delta-v | 1.184 m/s |
| Closed-loop dispersion runs | 5,000 |

These values describe the committed analysis outputs and their modeled assumptions. They are not flight-performance guarantees.

## Engineering approach

The project is organized around a central question: **Which mission constraints prevent a usable docking window, and what changes when those constraints are varied?**

### 1. Orbital and attitude modeling

The analysis establishes the target-orbit geometry, transforms states between relevant reference frames, and reconstructs the supplied attitude profile. It evaluates the effects of time handling, Earth orientation, quaternion conventions, and alternative attitude interpolation methods.

### 2. Constraint evaluation and event finding

The window-search pipeline evaluates the relevant conditions independently:

- Ground-station visibility and minimum elevation.
- Earth occultation and eclipse exclusion.
- Solar-port-angle constraints.
- Continuous-window duration.

Constraint boundaries are found through event detection and root refinement. The resulting feasible intervals are intersected to determine whether a candidate satisfies every requirement for the full specified duration.

### 3. Infeasibility diagnosis

Instead of returning only a binary feasible/infeasible answer, the analysis ranks near-miss candidates and records the constraints limiting each opportunity.

Sensitivity studies examine orbital phase, local time of the ascending node (LTAN), and alternative interpretations of underspecified mission parameters. The relaxation analysis then quantifies how far the illumination constraint must change before a candidate window reaches the required duration.

### 4. Conditional terminal-approach design

For the candidate window obtained under the tested relaxation, the project evaluates approach directions and constructs a nominal +V-bar approach using Clohessy–Wiltshire relative-motion dynamics.

The analysis includes maneuver targeting, hold points, final translation, missed-burn and passive-safety considerations, correction opportunities, and dispersion studies. A separate nonlinear J2 propagation provides a cross-check of the nominal relative trajectory.

### 5. Independent numerical checks

The validation pipeline compares selected calculations against independent astrodynamics references, including Astropy/ERFA, JPL DE421, and Orekit 13. It also evaluates event-finder grid sensitivity and records numerical comparisons in generated validation outputs.

## Recorded validation results

The existing run artifacts report the following checks:

| Check | Recorded result |
|---|---|
| Attitude-law identification | Maximum residual of approximately \(2.2 \times 10^{-11}\) rad across 14,401 samples |
| Sun-direction comparison | Approximately 0.01 arcsec difference against JPL DE421 |
| Reference-frame and station geometry | Approximately 2 cm position difference and \(5 \times 10^{-7}\)° elevation difference in the reported comparison |
| Event-finder grid sensitivity | Agreement within 0.7 ms across 1 s, 2 s, and 5 s grids |
| Independent end-to-end window check | Longest feasible candidate span agrees within 0.71 s |
| Orekit 13 comparison | Approximately 1.9 cm position difference over 20 days; AOS/LOS within 0.91 ms |
| Nonlinear J2 trajectory comparison | Relative deviation of approximately 5.1 cm or less at evaluated waypoints and contact |
| Built-in self-test | Nine checks passed in the recorded run |

These figures are specific to the recorded configuration and comparison procedures. Independent numerical agreement does not establish that every physical assumption is correct, and the results should be regenerated after material changes to the code or configuration.

## Repository structure

```text
aule-rpod-submission/
├── config/
│   └── mission.yaml
├── data/
├── outputs/
│   ├── part1/
│   ├── part2/
│   ├── sensitivity/
│   └── validation/
├── scripts/
├── src/
│   └── aule_rpod/
│       ├── orbit/
│       ├── attitude/
│       ├── geometry/
│       ├── events/
│       ├── rpod/
│       ├── validation/
│       ├── viz/
│       ├── windows.py
│       ├── sensitivity.py
│       ├── studies.py
│       ├── scenario.py
│       └── cli.py
├── tests/
├── Makefile
├── pyproject.toml
└── README.md
```

The tree above describes the principal components; consult the repository for the exact files included in the current revision.

## Reproducibility

### Requirements

- Python 3.10 or newer.
- Git.
- A Windows, Linux, or macOS environment.
- Additional validation data and dependencies for the full independent-validation workflow.

### Windows PowerShell

Run the following commands from the repository root:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Run the test suite and built-in self-test:

```powershell
python -m pytest
aule-rpod selftest
```

If PowerShell blocks virtual-environment activation, allow it for the current terminal session:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

### Linux and macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
aule-rpod selftest
```

### Run the analysis stages

```powershell
aule-rpod part1
aule-rpod part2
aule-rpod studies --quick
```

- `part1` runs the docking-window analysis.
- `part2` runs the conditional terminal-approach analysis.
- `studies --quick` runs the quick sensitivity and interpretation study.

The full analysis and independent-validation workflow may take longer and may require external reference data. Inspect the available command-line help and `Makefile` before running the full pipeline.

## Generated outputs

The analysis writes machine-readable tables, summaries, figures, and validation artifacts under `outputs/`.

| Output area | Purpose |
|---|---|
| `outputs/part1/` | Docking-window candidates, constraint events, and baseline feasibility results |
| `outputs/part2/` | Approach strategy, maneuver plan, trajectory analysis, and dispersion results |
| `outputs/sensitivity/` | Feasibility sweeps and sensitivity figures |
| `outputs/validation/` | Numerical comparison results and validation report |
| `outputs/run_manifest.json` | Recorded run provenance and configuration information |

The generated outputs should be interpreted together with the configuration used to produce them. A result from one mission scenario should not be assumed to generalize to another.

## Modeling assumptions and limitations

- The baseline uses a 550 km circular sun-synchronous orbit with secular J2 modeling; nonlinear J2 propagation is used for selected cross-checks.
- The analysis uses explicit reference-frame, time-scale, attitude-interpolation, and station-geometry conventions.
- The baseline constraints include a 30°–60° solar-port-angle band, eclipse exclusion, a minimum station elevation of 5°, and a continuous docking-window requirement of 540 seconds.
- The terminal-approach design uses Clohessy–Wiltshire targeting, with nonlinear propagation and dispersion studies used as additional checks.
- Monte Carlo results depend on the assumed initial-state errors, maneuver uncertainties, controller behavior, and other modeled parameters.
- The results do not establish flight readiness, flight-qualified guidance, hardware compatibility, or operational safety.
- Further work would require mission-specific requirements verification, more complete force and disturbance models where applicable, validated navigation and control assumptions, and independent review of all safety-critical constraints.

## Project objective

The goal is not merely to produce a trajectory or a plot. It is to make the feasibility decision auditable:

1. State the constraints.
2. Identify their event boundaries.
3. Determine whether a continuous feasible interval exists.
4. Explain why the best candidate fails when no interval qualifies.
5. Quantify the effect of changing an assumption.
6. Design and independently check a conditional approach when a candidate becomes available.

That separation between feasibility, sensitivity, and approach design is the core engineering structure of this project.

---

**Scope:** Research and mission-analysis prototype. Not certified for operational spacecraft use.
