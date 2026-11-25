"""Give each PASEF spectrum a unique Comet scan key while retaining its native ID."""

import argparse
import csv
import os
import re
from pathlib import Path


INDEX = re.compile(rb"\bindex=(\d+)\b")


def normalize(source, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    crosswalk = output.with_suffix(".scan_map.tsv")
    pending_output = output.with_suffix(".mgf.part")
    pending_map = crosswalk.with_suffix(".tsv.part")
    indexes, frames = set(), set()
    with source.open("rb") as data, pending_output.open("wb") as mgf, pending_map.open("w", newline="") as tsv:
        writer = csv.writer(tsv, delimiter="\t")
        writer.writerow(["comet_scan", "author_index", "native_scans", "rt_seconds", "precursor_mz", "charge"])
        block = []
        for line in data:
            if line.startswith(b"BEGIN IONS"):
                if block:
                    raise ValueError("Nested spectrum in input")
                block.append(line)
            elif block:
                block.append(line)
                if not line.startswith(b"END IONS"):
                    continue
                tags = {
                    key: (i, content.strip()) for i, row in enumerate(block)
                    if b"=" in row for key, content in [row.split(b"=", 1)]
                }
                match = INDEX.search(tags[b"TITLE"][1])
                if not match:
                    raise ValueError("MGF TITLE has no index")
                index = int(match.group(1))
                if index in indexes:
                    raise ValueError(f"Duplicate original spectrum index={index}")
                indexes.add(index)
                scan_line, native = tags[b"SCANS"]
                scan = index + 1
                block[scan_line] = f"SCANS={scan}\n".encode()
                mgf.writelines(block)
                frames.add(native)
                writer.writerow([
                    scan, index, native.decode(), tags[b"RTINSECONDS"][1].decode(),
                    tags[b"PEPMASS"][1].decode().split()[0], tags[b"CHARGE"][1].decode(),
                ])
                block.clear()
            elif line.strip():
                raise ValueError("Unexpected content outside an MGF spectrum")
        if block or not indexes:
            raise ValueError("Truncated or empty MGF input")
    os.replace(pending_output, output)
    os.replace(pending_map, crosswalk)
    print(f"{source.name}: {len(indexes)} spectra, {len(frames)} distinct native SCANS values", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    normalize(args.source, args.output)


if __name__ == "__main__":
    main()
