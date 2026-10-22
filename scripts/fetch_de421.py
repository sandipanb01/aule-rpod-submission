"""Fetch JPL DE421 Chebyshev arrays (PyPI 'de421' sdist; JPL/NASA public data) into data/de421/.

Used only by the optional V&V check of the Sun ephemeris. The 2008 package does not
build on modern Python, so the source archive is downloaded and its .npy arrays extracted.
"""
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "de421"


def main():
    if (DEST / "constants.npy").exists():
        print(f"DE421 already present in {DEST}")
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.check_call([sys.executable, "-m", "pip", "download", "--no-deps", "--no-binary", ":all:",
                               "-d", tmp, "de421==2008.1"])
        arc = next(Path(tmp).glob("de421-*.tar.gz"))
        DEST.mkdir(parents=True, exist_ok=True)
        with tarfile.open(arc) as tf:
            for m in tf.getmembers():
                if m.name.endswith(".npy"):
                    m.name = Path(m.name).name
                    tf.extract(m, DEST, filter="data")
    print(f"DE421 arrays extracted to {DEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
