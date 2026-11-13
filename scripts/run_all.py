"""Run the whole analysis pipeline in order (about 3-4 minutes).

    uv run scripts/run_all.py            # run steps 01-11
    uv run scripts/run_all.py --fresh    # delete results/ and figures/ first
    uv run scripts/run_all.py 05 09      # run only these steps (they read earlier results from results/)

Data must already be downloaded (scripts/download_data.py).
"""
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
steps = sorted(p for p in (ROOT / "scripts").glob("[0-9][0-9]_*.py"))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
if "--fresh" in sys.argv:
    for d in ("results", "figures"):
        shutil.rmtree(ROOT / d, ignore_errors=True)
if args:
    steps = [p for p in steps if p.name[:2] in args]
t0 = time.time()
for p in steps:
    t = time.time()
    r = subprocess.run([sys.executable, "-u", str(p)], cwd=ROOT)
    if r.returncode:
        sys.exit(f"step {p.name} failed")
    print(f"   ({time.time() - t:.0f} s)\n")
print(f"done: {len(steps)} steps in {time.time() - t0:.0f} s")
