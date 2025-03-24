#!/bin/bash

# Convert bedpe file as arg 1 into fragments file as arg 1
set -euo pipefail

SAMPLE_ID=$1
TMP_DIR="$2/$1"
OUT_PATH=$3

mkdir -p $TMP_DIR
URL="http://catlas.org/catlas_downloads/humanbrain/bedpe/$1.bedpe.gz"

gcc brain/count_unique.c -O3 -o $TMP_DIR/count_unique

curl "$URL" |
    gunzip -c | \
    python3 brain/convert_fragments.py |
    LC_ALL=C sort --parallel=4 -k1,1Vf -k2,2n -S 1G |
    $TMP_DIR/count_unique |
    bgzip -@ 4 > $OUT_PATH;
tabix -p bed $OUT_PATH


