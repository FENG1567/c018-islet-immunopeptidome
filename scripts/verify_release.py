#!/usr/bin/env python3
"""Check the hashes of the files distributed in this release."""

import argparse
import csv
import hashlib
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.tsv"
FIELDS = ("path", "bytes", "sha256")


def digest(file):
    with file.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_manifest():
    excluded = {".git", ".venv", "__pycache__", "build", "input"}
    with MANIFEST.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for file in sorted(ROOT.rglob("*")):
            relative = file.relative_to(ROOT)
            if file == MANIFEST or excluded.intersection(relative.parts) or not file.is_file():
                continue
            if file.is_symlink():
                raise ValueError(f"Release files must be regular files: {relative}")
            writer.writerow(
                dict(path=relative.as_posix(), bytes=file.stat().st_size, sha256=digest(file))
            )


def validate_manifest():
    seen = set()
    with MANIFEST.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Unexpected release manifest header")
        for row in reader:
            relative = PurePosixPath(row["path"])
            if relative.is_absolute() or ".." in relative.parts or relative in seen:
                raise ValueError(f"Unsafe or duplicate release path: {relative}")
            seen.add(relative)
            file = ROOT.joinpath(*relative.parts)
            if (not file.is_file() or not file.resolve().is_relative_to(ROOT)
                    or file.stat().st_size != int(row["bytes"]) or digest(file) != row["sha256"]):
                raise ValueError(f"Release file is missing or has changed: {relative}")
    if not seen:
        raise ValueError("Empty release manifest")
    print(f"Verified {len(seen)} release files.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true", help="Freeze a release after authoring.")
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest()
    validate_manifest()


if __name__ == "__main__":
    main()
