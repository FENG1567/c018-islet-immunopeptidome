"""Apply the fixed matched-panel search as verified source files become available."""

import argparse
import subprocess
import sys
import time
from pathlib import Path

from matched_panel_search import search_one


RUN_SIZES = {
    "20180420_Stern_Islets_W632": 13185889,
    "20201202_Lumos_Stern_Islet_Set1_W632_5": 29769383,
    "20201202_Lumos_Stern_Islet_Set2_W632_10": 22595113,
    "20200228_IsletsW632_OTOT_lumos": 44988769,
    "20200302_IsletsW6322_OTOT_lumos": 41413600,
    "20200724_Lumos_Stern_W632_OTOT_5": 28246794,
    "20190405_Stern_Islet_Control_W632": 35322367,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--comet", type=Path, required=True, help="Comet executable")
    parser.add_argument("--watch-seconds", type=int, default=0)
    args = parser.parse_args()
    base, source, output = (p.resolve(strict=True) for p in (args.base, args.source, args.output))
    deadline = time.monotonic() + args.watch_seconds
    completed = set()
    while True:
        for run, size in RUN_SIZES.items():
            mgf = source / (run + ".mgf")
            if run not in completed and mgf.is_file() and mgf.stat().st_size == size:
                search_one(base, output, mgf, args.comet.resolve(strict=True))
                completed.add(run)
        if len(completed) == len(RUN_SIZES) or time.monotonic() >= deadline:
            break
        time.sleep(20)
    if completed:
        pepxml = [str(output / "body" / (run + ".pep.xml")) for run in sorted(completed)]
        subprocess.run([sys.executable, "-u", str(Path(__file__).with_name("score_panel.py")),
                        "--panel", str(output / "reference_targets_shuffles.tsv"), "--out",
                        str(output / "summary"), "--pepxml", *pepxml], check=True)
    print("matched panel scored", sorted(completed), flush=True)


if __name__ == "__main__":
    main()
