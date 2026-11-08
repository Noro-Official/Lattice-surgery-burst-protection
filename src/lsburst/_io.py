"""File, checksum and comparison helpers shared by the analyses and tools."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"


def sha256_file(path: Path) -> str:
    """Returns the SHA-256 hex digest of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_checksum_list(path: Path) -> dict[str, str]:
    """Reads a `sha256sum`-style list into {relative path: digest}.

    Accepts both `digest  path` and `digest *./path` line formats.
    """
    sums = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(maxsplit=1)
        sums[name.strip().lstrip("*").removeprefix("./")] = digest
    return sums


def load_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def to_jsonable(obj: Any) -> Any:
    """Converts numpy scalars and arrays inside `obj` to plain Python values.

    Round-trips through `json` so that the result is exactly what would be read back from disk.
    """
    def default(x):
        return x.tolist() if isinstance(x, np.ndarray) else float(x)
    return json.loads(json.dumps(obj, default=default))


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(to_jsonable(obj), indent=1))


def tree_differences(expected: Any, actual: Any, path: str = "", rel_tol: float = 1e-9,
                     skip_keys: tuple[str, ...] = ("provenance", "authenticity", "source")) -> list[str]:
    """Lists the places where two JSON-like trees differ.

    Every key of `expected` must be present in `actual`. Floats are compared with a relative tolerance of
    `rel_tol` (absolute 1e-15); NaN equals NaN. Keys in `skip_keys` are ignored.

    Args:
        expected: The stored (reference) tree.
        actual: The recomputed tree.
        path: Prefix used in the returned messages.
        rel_tol: Relative tolerance for floats.
        skip_keys: Metadata keys that are not part of the analysis.

    Returns:
        One message per difference; empty when the trees agree.
    """
    diffs: list[str] = []
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return [f"{path}: expected an object"]
        for key, value in expected.items():
            if key in skip_keys:
                continue
            if key not in actual:
                diffs.append(f"{path}/{key}: missing")
            else:
                diffs += tree_differences(value, actual[key], f"{path}/{key}", rel_tol, skip_keys)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            return [f"{path}: list length differs"]
        for i, (x, y) in enumerate(zip(expected, actual)):
            diffs += tree_differences(x, y, f"{path}[{i}]", rel_tol, skip_keys)
    elif isinstance(expected, float) and isinstance(actual, (int, float)):
        same = math.isclose(expected, actual, rel_tol=rel_tol, abs_tol=1e-15)
        if not (same or (math.isnan(expected) and math.isnan(actual))):
            diffs.append(f"{path}: {expected!r} != {actual!r}")
    elif expected != actual:
        diffs.append(f"{path}: {expected!r} != {actual!r}")
    return diffs
