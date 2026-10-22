"""Reader/writer for STK attitude files (*.a), AttitudeTimeQuaternions format.

STK convention (Ansys STK Help, 'Attitude File Format (*.a)'):
  * "The quaternions represent the rotation from the reference coordinate
    frame to the vehicle's body coordinate frame."
  * "The fourth element of the quaternion is the scalar component."
  * CoordinateAxes = reference axes of the rotation (here: LVLH of the target).
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class AttitudeProfile:
    times_s: np.ndarray            # seconds from ScenarioEpoch
    quat_xyzw: np.ndarray          # (N, 4) scalar-last
    scenario_epoch: str
    coordinate_axes: str
    meta: dict = field(default_factory=dict)

    @property
    def step_s(self) -> float:
        d = np.diff(self.times_s)
        return float(np.median(d))


_MONTHS = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}


def stk_date_to_iso(s: str) -> str:
    """'1 Jul 2026 06:30:00.000000' -> '2026-07-01T06:30:00.000000'."""
    d, mon, y, hms = s.split()[:4]
    return f"{int(y):04d}-{_MONTHS[mon[:3]]:02d}-{int(d):02d}T{hms}"


def read_stk_attitude(path: str | Path) -> AttitudeProfile:
    path = Path(path)
    raw = path.read_bytes()
    meta = {"file": path.name, "sha256": hashlib.sha256(raw).hexdigest()}
    epoch, axes, n_decl = None, None, None
    rows = []
    in_block = False
    for line in raw.decode("ascii", errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            if s.startswith("# WrittenBy"):
                meta["written_by"] = s.split()[-1]
            continue
        if s.upper().startswith("BEGIN ATTITUDE"):
            in_block = True
            continue
        if s.upper().startswith("END ATTITUDE"):
            break
        if not in_block:
            meta.setdefault("version", s)
            continue
        key = s.split()[0]
        if key.lower() == "scenarioepoch":
            epoch = stk_date_to_iso(s.split(None, 1)[1])
        elif key.lower() == "coordinateaxes":
            axes = " ".join(s.split()[1:])
        elif key.lower() == "numberofattitudepoints":
            n_decl = int(s.split()[1])
        elif key.lower() in ("interpolationmethod", "interpolationorder", "centralbody",
                             "sequence", "blockingfactor"):
            meta[key] = s.split()[1]
        else:
            toks = s.split()
            if len(toks) == 5 and all(re.match(r"^[-+0-9.eE]+$", t) for t in toks):
                rows.append([float(t) for t in toks])
    if not rows:
        raise ValueError(f"No attitude records found in {path}")
    arr = np.asarray(rows)
    q = arr[:, 1:5]
    norm = np.linalg.norm(q, axis=1)
    meta.update({"declared_points": n_decl, "parsed_points": len(arr),
                 "max_norm_error": float(np.max(np.abs(norm - 1.0)))})
    if n_decl is not None and n_decl != len(arr):
        raise ValueError(f"Declared {n_decl} points but parsed {len(arr)}")
    q = q / norm[:, None]                       # STK normalises on read
    # enforce hemisphere continuity (q and -q are the same attitude)
    for k in range(1, len(q)):
        if np.dot(q[k], q[k - 1]) < 0:
            q[k] = -q[k]
    return AttitudeProfile(arr[:, 0], q, epoch, axes, meta)


def write_stk_attitude(path, times_s, quat_xyzw, scenario_epoch_stk: str, axes: str = "LVLH Satellite/Target"):
    with open(path, "w") as fh:
        fh.write("stk.v.12.0\n\nBEGIN Attitude\n")
        fh.write(f"    NumberOfAttitudePoints {len(times_s)}\n")
        fh.write(f"    ScenarioEpoch {scenario_epoch_stk}\n")
        fh.write(f"    CoordinateAxes {axes}\n")
        fh.write("    AttitudeTimeQuaternions\n")
        for t, q in zip(times_s, quat_xyzw):
            fh.write(f" {t:.16e} {q[0]: .16e} {q[1]: .16e} {q[2]: .16e} {q[3]: .16e}\n")
        fh.write("END Attitude\n")
