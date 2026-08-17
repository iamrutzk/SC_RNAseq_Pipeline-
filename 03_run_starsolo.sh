#!/usr/bin/env bash
# Phase 3.1 (terminal version of "RNA STARsolo" Galaxy tool).
# Usage: 03_run_starsolo.sh config/config.yaml <sample_name> <R1_cDNA.fastq.gz,...> <R2_barcode.fastq.gz,...>
#
# NOTE on read order: STARsolo wants cDNA reads first, then barcode reads:
#   --readFilesIn cDNA_R2_lanes_comma_separated  barcode_R1_lanes_comma_separated
# This mirrors the Galaxy tool's "Barcode reads" vs "cDNA reads" split.
set -euo pipefail

CONFIG="${1:?config path}"
SAMPLE="${2:?sample name}"
CDNA_READS="${3:?comma-separated cDNA fastq.gz files}"
BARCODE_READS="${4:?comma-separated barcode fastq.gz files}"

REF_DIR=$(yq -r '.paths.ref_dir' "$CONFIG")
RESULTS_DIR=$(yq -r '.paths.results_dir' "$CONFIG")
CHEMISTRY=$(yq -r '.chemistry' "$CONFIG")
THREADS=$(yq -r '.reference.star_index_threads' "$CONFIG")

OUTDIR="${RESULTS_DIR}/${SAMPLE}/starsolo"
mkdir -p "$OUTDIR"

case "$CHEMISTRY" in
  10x_v3)
    CB_LEN=$(yq -r '.solo_10x_v3.cb_len' "$CONFIG")
    UMI_LEN=$(yq -r '.solo_10x_v3.umi_len' "$CONFIG")
    WL_URL=$(yq -r '.solo_10x_v3.whitelist_url' "$CONFIG")
    WHITELIST="${REF_DIR}/$(basename "${WL_URL%.gz}")"
    ;;
  10x_v2)
    CB_LEN=$(yq -r '.solo_10x_v2.cb_len' "$CONFIG")
    UMI_LEN=$(yq -r '.solo_10x_v2.umi_len' "$CONFIG")
    WL_URL=$(yq -r '.solo_10x_v2.whitelist_url' "$CONFIG")
    WHITELIST="${REF_DIR}/$(basename "${WL_URL%.gz}")"
    ;;
  indrop)
    echo "!! chemistry=indrop: verified CB/UMI start,len values are NOT set" >&2
    echo "!! by default. Run scripts/02b_inspect_read_structure.sh first," >&2
    echo "!! fill in CB_LEN/UMI_LEN/WHITELIST below, then re-run." >&2
    exit 1
    ;;
  *)
    echo "Unknown chemistry '$CHEMISTRY' in config.yaml" >&2
    exit 1
    ;;
esac

echo ">>> Running STARsolo for sample: $SAMPLE (chemistry=$CHEMISTRY)"

STAR \
  --runMode alignReads \
  --genomeDir "${REF_DIR}/star_index" \
  --readFilesIn "$CDNA_READS" "$BARCODE_READS" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WHITELIST" \
  --soloCBstart 1 --soloCBlen "$CB_LEN" \
  --soloUMIstart $((CB_LEN + 1)) --soloUMIlen "$UMI_LEN" \
  --soloCBmatchWLtype 1MM_multi \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloFeatures Gene \
  --soloCellFilter None \
  --outSAMattributes NH HI AS nM CB UB \
  --outSAMtype BAM Unsorted \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "${OUTDIR}/"

echo ">>> STARsolo finished. Raw matrix at:"
echo "    ${OUTDIR}/Solo.out/Gene/raw/{matrix.mtx,barcodes.tsv,features.tsv}"
echo ">>> Mapping QC summary:"
cat "${OUTDIR}/Solo.out/Gene/Summary.csv" || true
