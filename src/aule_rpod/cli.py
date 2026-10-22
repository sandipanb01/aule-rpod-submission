"""Single command-line entry point.

    aule-rpod selftest            built-in health check (exit code 1 on failure)
    aule-rpod part1               docking-window analysis -> outputs/part1
    aule-rpod part2               approach design at the first window -> outputs/part2
    aule-rpod studies [--quick]   sensitivity / alternatives -> outputs/sensitivity
    aule-rpod validate            independent V&V (astropy, DE421, Orekit) -> outputs/validation
    aule-rpod all                 everything above, plus a provenance manifest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from importlib import metadata
from pathlib import Path

from .config import REPO_ROOT, load_config


def _manifest(cfg, outdir: Path, steps):
    src = sorted((REPO_ROOT / "src" / "aule_rpod").rglob("*.py"))
    h = hashlib.sha256()
    for f in src:
        h.update(f.read_bytes())
    vers = {}
    for pkg in ("numpy", "scipy", "pyerfa", "astropy", "astropy-iers-data", "pandas", "matplotlib", "plotly",
                "pydantic", "orekit-jpype"):
        try:
            vers[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            vers[pkg] = None
    att = cfg.attitude_path.read_bytes()
    m = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "python": platform.python_version(),
         "platform": platform.platform(), "config_fingerprint": cfg.fingerprint(),
         "source_sha256": h.hexdigest(), "source_files": len(src),
         "input_attitude_sha256": hashlib.sha256(att).hexdigest(), "packages": vers, "steps": steps}
    (outdir / "run_manifest.json").write_text(json.dumps(m, indent=2))
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(prog="aule-rpod", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "part1", "part2", "studies", "validate", "all"])
    ap.add_argument("--config", default=None, help="mission YAML (default config/mission.yaml)")
    ap.add_argument("--out", default=str(REPO_ROOT / "outputs"), help="output root")
    ap.add_argument("--quick", action="store_true", help="coarse grids for the sensitivity maps")
    ap.add_argument("--de421", default=str(REPO_ROOT / "data" / "de421"),
                    help="directory of JPL DE421 arrays (scripts/fetch_de421.py; optional V&V)")
    ap.add_argument("--orekit-data", default=str(REPO_ROOT / "data" / "orekit-data"))
    ap.add_argument("--mc-samples", type=int, default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    steps = {}
    t0 = time.time()
    if a.command in ("selftest", "all"):
        from . import selftest
        res, _ = selftest.run(cfg)
        steps["selftest"] = res.value
        if res.value != "PASSED":
            return 1
    if a.command in ("part1", "all"):
        from . import part1
        r = part1.run(cfg, out / "part1")
        steps["part1"] = {"strict_windows": r["summary"]["strict_docking_windows"],
                          "first_window_min_relaxation": r["summary"]["first_window_under_min_relaxation_utc"]}
    if a.command in ("part2", "all"):
        from . import part2
        r = part2.run(cfg, out / "part2", mc_samples=a.mc_samples)
        steps["part2"] = {"strategy": r["summary"]["selected_strategy"], "dv_total_mps": r["summary"]["dv_total_mps"]}
    if a.command in ("studies", "all"):
        from . import studies
        studies.run(cfg, out / "sensitivity", quick=a.quick)
        steps["studies"] = "done"
    if a.command in ("validate", "all"):
        from .validation import run as V
        cand = out / "part1" / "csv" / "docking_window_candidates.csv"
        rep = V.run(cfg, out / "validation", de421_dir=a.de421, orekit_data=a.orekit_data,
                    candidates_csv=cand if cand.exists() else None)
        steps["validate"] = {k: v.get("verdict") for k, v in rep.items() if isinstance(v, dict)}
    steps["wall_time_s"] = round(time.time() - t0, 1)
    if a.command == "all":
        _manifest(cfg, out, steps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
