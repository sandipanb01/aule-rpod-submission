"""Build a minimal Orekit data directory (data/orekit-data) from bundled sources:
  * tai-utc.dat  (USNO format) from ERFA's IERS leap-second table
  * finals2000A.all (IERS EOP, measured values) from the astropy-iers-data package
Orekit is only used for independent V&V; nothing here touches the main pipeline.
"""
import shutil
import sys
from pathlib import Path

import erfa

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "orekit-data"
MON = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    lines = []
    for row in erfa.leap_seconds.get():
        y, m, dat = int(row["year"]), int(row["month"]), float(row["tai_utc"])
        if y < 1972:
            continue
        jd = sum(erfa.cal2jd(y, m, 1))
        lines.append(f" {y:4d} {MON[m - 1]:>3s}  1 =JD {jd:9.1f}  TAI-UTC=  {dat:10.7f} S + (MJD - 41317.) X 0.0      S")
    (DEST / "tai-utc.dat").write_text("\n".join(lines) + "\n")
    import astropy_iers_data
    src = Path(astropy_iers_data.__file__).parent / "data" / "finals2000A.all"
    shutil.copy(src, DEST / "finals2000A.all")
    print(f"Orekit data written to {DEST} ({len(lines)} leap-second records + finals2000A.all)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
