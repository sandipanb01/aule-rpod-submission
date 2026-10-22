"""Assemble the submission archive: source, configuration, tests, input data, the
generated products and the typeset report (Orekit data and DE421 arrays are left
out; `make data` rebuilds them from public sources).

    python scripts/package_submission.py [--out DIR]
"""
from __future__ import annotations

import argparse
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INCLUDE = ["README.md", "Makefile", "pyproject.toml", "requirements.lock", "config", "src", "tests", "scripts",
           "data/Target_Attitude_Profile.a", "docs/report/report.pdf", "docs/report/report_source.md",
           "docs/report/build_pdf.py", "outputs"]
SKIP_DIRS = {"__pycache__", ".pytest_cache", "build", ".ipynb_checkpoints"}
SKIP_SUFFIX = {".pyc", ".egg-info"}


def files():
    for item in INCLUDE:
        p = ROOT / item
        if p.is_file():
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.is_file() and not (set(f.relative_to(ROOT).parts) & SKIP_DIRS) \
                        and f.suffix not in SKIP_SUFFIX and not any(part.endswith(".egg-info") for part in f.parts):
                    yield f
        else:
            raise FileNotFoundError(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "dist"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    name = "aule-rpod-submission"
    zpath = out / f"{name}.zip"
    n = 0
    with zipfile.ZipFile(zpath, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for f in files():
            z.write(f, f"{name}/{f.relative_to(ROOT).as_posix()}")
            n += 1
        log = ROOT / "outputs_full_run.log"
        if log.exists():
            z.write(log, f"{name}/outputs/run_log.txt")
            n += 1
    shutil.copy2(ROOT / "docs/report/report.pdf", out / "Aule_RPOD_Technical_Report.pdf")
    print(f"{zpath} ({n} files, {zpath.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
