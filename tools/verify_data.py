#!/usr/bin/env python3
"""Verifies the SHA-256 checksums of every data file.

Two lists are checked:
    SHA256SUMS.txt                           every file under data/, written for this repository
    data/p1_archive/PACKAGE_SHA256SUMS.txt   the 31 October 2026 publication package's own checksums for the archived
                                             cohort files (per-event blocks are checked when they are loaded, in
                                             step 2)

Usage:
    python tools/verify_data.py
"""

from __future__ import annotations

import argparse
import pathlib
import sys

src_path = pathlib.Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(src_path))

from lsburst._io import DATA_DIR, REPO_ROOT, read_checksum_list, sha256_file  # noqa: E402

PACKAGE_PREFIX = "experiments/pacer_p1/results/"


def check_repository_list(path: pathlib.Path) -> tuple[int, list[str]]:
    """Returns (files checked, files that are missing or differ) for SHA256SUMS.txt."""
    sums = read_checksum_list(path)
    bad = [name for name, digest in sums.items()
           if not (REPO_ROOT / name).is_file() or sha256_file(REPO_ROOT / name) != digest]
    listed = set(sums)
    unlisted = [str(p.relative_to(REPO_ROOT)) for p in DATA_DIR.rglob("*")
                if p.is_file() and str(p.relative_to(REPO_ROOT)) not in listed]
    return len(sums), bad + [f"{name} (not listed)" for name in unlisted]


def check_package_list(archive_dir: pathlib.Path) -> tuple[int, list[str]]:
    """Checks the archived cohort files (not the bundled blocks) against the package's own checksums."""
    sums = read_checksum_list(archive_dir / "PACKAGE_SHA256SUMS.txt")
    files = [p for p in archive_dir.rglob("*")
             if p.is_file() and p.suffix != ".md" and p.name not in ("PACKAGE_SHA256SUMS.txt", "blocks.jsonl.gz")]
    bad = [str(p.relative_to(archive_dir)) for p in files
           if sums.get(PACKAGE_PREFIX + str(p.relative_to(archive_dir))) != sha256_file(p)]
    return len(files), bad


def main() -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    n_repo, bad_repo = check_repository_list(REPO_ROOT / "SHA256SUMS.txt")
    n_pkg, bad_pkg = check_package_list(DATA_DIR / "p1_archive")
    print(f"SHA256SUMS.txt: {n_repo - len(bad_repo)}/{n_repo} files verified")
    print(f"PACKAGE_SHA256SUMS.txt: {n_pkg - len(bad_pkg)}/{n_pkg} archived cohort files match the original package")
    for name in bad_repo + bad_pkg:
        print(f"  MISMATCH {name}", file=sys.stderr)
    return 1 if bad_repo or bad_pkg else 0


if __name__ == "__main__":
    sys.exit(main())
