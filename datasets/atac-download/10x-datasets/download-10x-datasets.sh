set -euo pipefail



OUTPUT_DIR=$1/atac

mkdir -p $OUTPUT_DIR


mkdir -p $OUTPUT_DIR/3k_pbmc/fragments
curl -Lo $OUTPUT_DIR/3k_pbmc/fragments/3k_pbmc.fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
curl -Lo $OUTPUT_DIR/3k_pbmc/fragments/3k_pbmc.fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi
curl -Lo $OUTPUT_DIR/3k_pbmc/peaks.bed https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_peaks.bed
sed -i '/^#/d' $OUTPUT_DIR/3k_pbmc/peaks.bed # Strip comment lines
Rscript fragments-to-bpcells.R $OUTPUT_DIR/3k_pbmc/fragments/3k_pbmc.fragments.tsv.gz $OUTPUT_DIR/3k_pbmc/bpcells/3k_pbmc



mkdir -p $OUTPUT_DIR/10k_pbmc/fragments
curl -Lo $OUTPUT_DIR/10k_pbmc/fragments/10k_pbmc.fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz
curl -Lo $OUTPUT_DIR/10k_pbmc/fragments/10k_pbmc.fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi
curl -Lo $OUTPUT_DIR/10k_pbmc/peaks.bed https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_peaks.bed
sed -i '/^#/d' $OUTPUT_DIR/10k_pbmc/peaks.bed # Strip comment lines
Rscript fragments-to-bpcells.R $OUTPUT_DIR/10k_pbmc/fragments/10k_pbmc.fragments.tsv.gz $OUTPUT_DIR/10k_pbmc/bpcells/10k_pbmc