set -euo pipefail

OUT_DIR=$1

mkdir -p $OUT_DIR/45k_pancreas
mkdir -p $OUT_DIR/500k_heart

# Download peak sets
curl -L https://www.encodeproject.org/files/ENCFF678PAI/@@download/ENCFF678PAI.bed.gz | \
    gunzip -c | \
    cut -f 1-3 | \
    sort -k1,1V -k2,2n -u > \
    $OUT_DIR/45k_pancreas/peaks.bed

curl -L https://www.encodeproject.org/files/ENCFF012BKZ/@@download/ENCFF012BKZ.bed.gz | \
    gunzip -c | \
    cut -f 1-3 | \
    sort -k1,1V -k2,2n -u > \
    $OUT_DIR/500k_heart/peaks.bed