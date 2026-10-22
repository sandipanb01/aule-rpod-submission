"""Typed, validated mission configuration (schema-first, units in field names).

Every product the pipeline writes is traceable to one YAML file. The schema
rejects unknown keys and out-of-range values at load time, so a typo cannot
silently change an analysis.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import List, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

REPO_ROOT = Path(__file__).resolve().parents[2]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class Scenario(_Strict):
    name: str
    orbit_epoch_utc: str
    mission_start_utc: str
    mission_end_utc: str


class TargetOrbit(_Strict):
    altitude_km: float = Field(gt=100.0, lt=5000.0)
    eccentricity: float = Field(ge=0.0, lt=0.1)
    ltan_hours: float = Field(ge=0.0, lt=24.0)
    arg_latitude_at_epoch_deg: float = 0.0
    propagator: Literal["secular_j2", "cowell_osculating", "cowell_mean"] = "secular_j2"


class Earth(_Strict):
    mu_km3_s2: float
    equatorial_radius_km: float
    flattening: float
    j2: float
    tropical_year_days: float


class Attitude(_Strict):
    file: str
    interpolation: Literal["slerp", "nlerp", "law"] = "slerp"
    docking_port_body: List[float]

    @field_validator("docking_port_body")
    @classmethod
    def _three(cls, v):
        if len(v) != 3:
            raise ValueError("docking_port_body must have 3 components")
        return v


class Constraints(_Strict):
    sun_angle_min_deg: float = Field(ge=0.0, le=180.0)
    sun_angle_max_deg: float = Field(ge=0.0, le=180.0)
    docking_duration_s: float = Field(gt=0.0)
    min_elevation_deg: float = Field(ge=-5.0, le=90.0)
    eclipse_excludes_penumbra: bool = True


class GroundStation(_Strict):
    name: str
    lat_deg: float = Field(ge=-90.0, le=90.0)
    lon_deg: float = Field(ge=-180.0, le=360.0)
    alt_m: float = 0.0


class Numerics(_Strict):
    grid_step_s: float = Field(gt=0.0, le=30.0)
    root_tol_s: float = Field(gt=0.0, le=1.0)
    timeseries_step_s: float = Field(gt=0.0)


class MonteCarlo(_Strict):
    samples: int = Field(ge=100)
    seed: int
    nav_position_sigma_frac_of_range: float
    nav_position_sigma_floor_m: float
    nav_velocity_sigma_mps: float
    burn_magnitude_sigma_frac: float
    burn_pointing_sigma_deg: float
    nav_correlation_time_s: float = Field(default=30.0, gt=0)
    control_period_s: float = Field(default=2.0, gt=0)
    lateral_time_constant_s: float = Field(default=20.0, gt=0)
    lateral_speed_limit_mps: float = Field(default=0.02, gt=0)
    thruster_deadband_mps: float = Field(default=0.002, ge=0)


class Part2(_Strict):
    initial_relative_position_lvlh_m: List[float]
    initial_relative_velocity_lvlh_mps: List[float]
    hold_point_distance_m: float
    hold_duration_s: float
    final_approach_speed_mps: float
    contact_speed_mps: float
    keep_out_sphere_radius_m: float
    collision_radius_m: float
    corridor_half_angle_deg: float
    camera_sun_exclusion_deg: float
    window_margin_s: float
    monte_carlo: MonteCarlo


class MissionConfig(_Strict):
    scenario: Scenario
    target_orbit: TargetOrbit
    earth: Earth
    attitude: Attitude
    constraints: Constraints
    ground_stations: List[GroundStation]
    numerics: Numerics
    part2: Part2

    @property
    def attitude_path(self) -> Path:
        p = Path(self.attitude.file)
        return p if p.is_absolute() else REPO_ROOT / p

    def fingerprint(self) -> str:
        """SHA-256 of the canonical JSON of the configuration (provenance)."""
        blob = json.dumps(self.model_dump(), sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:16]

    def with_updates(self, **sections) -> "MissionConfig":
        """Return a copy with selected sections/fields replaced, e.g.
        ``cfg.with_updates(target_orbit={'ltan_hours': 22.5})``."""
        data = self.model_dump()
        for sec, upd in sections.items():
            if isinstance(upd, dict):
                data[sec] = {**data[sec], **upd}
            else:
                data[sec] = upd
        return MissionConfig.model_validate(data)


def load_config(path: str | Path | None = None) -> MissionConfig:
    path = Path(path) if path else REPO_ROOT / "config" / "mission.yaml"
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return MissionConfig.model_validate(raw)
