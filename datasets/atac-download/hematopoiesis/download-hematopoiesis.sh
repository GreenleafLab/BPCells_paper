#!/bin/bash
#SBATCH --job-name=10x-download-atac
#SBATCH --partition=wjg,biochem,sfgf

set -euo pipefail

TMP_DIR=$1/35k_hematopoiesis
OUTPUT_DIR=$2
REFERENCE_DIR=$3

mkdir -p $TMP_DIR
mkdir -p $OUTPUT_DIR/fragments
mkdir -p $OUTPUT_DIR/cell-barcodes

# Download big tar of fragment files from GEO
if [ ! -f $TMP_DIR/hematopoiesis.tar ]; then
    curl -Lo $TMP_DIR/hematopoiesis.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139369&format=file&id=32%2C33%2C34%2C35%2C36%2C37%2C38%2C39%2C40%2C41'
fi

# Untar the data if the output dir is empty
mkdir -p $TMP_DIR/hematopoiesis
tar -xf $TMP_DIR/hematopoiesis.tar --dir $TMP_DIR/hematopoiesis

# Run tabix on all files in parallel (10 files)
# Also run data conversion + cell filtering
for f in $TMP_DIR/hematopoiesis/*.tsv.gz; do
    FILE=$(basename $f)
    FRAGMENT_NAME=${FILE#*_scATAC_}
    
    cp $f $OUTPUT_DIR/fragments/$FRAGMENT_NAME;
    tabix -p bed $OUTPUT_DIR/fragments/$FRAGMENT_NAME &
    set -x
    SAMPLE=${FRAGMENT_NAME%.fragments.tsv.gz}
    ( \
        Rscript fragments-to-bpcells.R $OUTPUT_DIR/fragments/$FRAGMENT_NAME $OUTPUT_DIR/bpcells/$SAMPLE \
        && Rscript find-passing-cells.R $OUTPUT_DIR/bpcells/$SAMPLE $OUTPUT_DIR/cell-barcodes/$SAMPLE.txt hg19 2000 10 $REFERENCE_DIR \
    ) &
    set +x
done
wait;

# Check all the files were created
ERROR=false
for f in $TMP_DIR/hematopoiesis/*.tsv.gz; do
    FILE=$(basename $f)
    FRAGMENT_NAME=${FILE#*_scATAC_}
    SAMPLE=${FRAGMENT_NAME%.fragments.tsv.gz}
    if [ ! -f $OUTPUT_DIR/fragments/$FRAGMENT_NAME.tbi ]; then
        echo "Missing file: $OUTPUT_DIR/fragments/$FRAGMENT_NAME.tbi"
        ERROR=true
    fi
    if [ ! -f $OUTPUT_DIR/cell-barcodes/$SAMPLE.txt ]; then
        echo "Missing file: $OUTPUT_DIR/cell-barcodes/$SAMPLE.txt"
        ERROR=true
    fi
done

if [ $ERROR == "true" ]; then
    exit 1
fi


# Download and convert peak sets
if [ ! -f $TMP_DIR/hematopoiesis_peaks.rds ]; then
    curl -Lo $TMP_DIR/hematopoiesis_peaks.rds 'https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Disease-Data/scATAC-All-Hematopoiesis-MPAL-191120.rds'
fi

Rscript hematopoiesis/convert_peaks.R  $TMP_DIR/hematopoiesis_peaks.rds $OUTPUT_DIR/peaks.bed

