#!/bin/bash
set -euo pipefail

DATASET_DIR=$1
TMP_DIR=$2
SAMPLE=$3

mkdir -p $TMP_DIR
cp "$DATASET_DIR/fragments/$SAMPLE.fragments.tsv.gz" "$TMP_DIR/sample.fragments.tsv.gz"
cp "$DATASET_DIR/fragments/$SAMPLE.fragments.tsv.gz.tbi" "$TMP_DIR/sample.fragments.tsv.gz.tbi"