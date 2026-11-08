#!/usr/bin/env bash
# Step 1: verify the SHA-256 checksums of every data file.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

python tools/verify_data.py
