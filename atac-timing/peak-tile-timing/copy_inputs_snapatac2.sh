#!/bin/bash
set -euo pipefail

DATASET_DIR="$1"
TMP_DIR="$2"
SAMPLE="$3"

mkdir -p $TMP_DIR
cp -r "$DATASET_DIR/snapatac2/$SAMPLE.filtered.h5ad" "$TMP_DIR/sample.h5ad"
