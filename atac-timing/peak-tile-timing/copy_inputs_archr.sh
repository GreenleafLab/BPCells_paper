#!/bin/bash
set -euo pipefail

DATASET_DIR=$1
TMP_DIR=$2
SAMPLE=$3

mkdir -p $TMP_DIR
cp "$DATASET_DIR/archr/$SAMPLE.filtered.arrow" "$TMP_DIR/sample.arrow"