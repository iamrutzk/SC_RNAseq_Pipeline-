#!/usr/bin/env bash
# Phase 2 (terminal version): pull FASTQ from NCBI SRA with sra-tools.
# Usage: 01_download_sra.sh config/config.yaml
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
SRA_DIR=$(yq -r '.paths.sra_dir' "$CONFIG")
FASTQ_DIR=$(yq -r '.paths.fastq_dir' "$CONFIG")
BIOPROJECT=$(yq -r '.bioproject' "$CONFIG")

mkdir -p "$SRA_DIR" "$FASTQ_DIR"

echo ">>> Verifying accession list for BioProject $BIOPROJECT"
echo "    (double-check data/sra/runinfo.csv against config.yaml's sra_accessions"
echo "     before trusting the placeholders shipped in this config)"

if command -v esearch >/dev/null 2>&1; then
  esearch -db sra -query "$BIOPROJECT" | efetch -format runinfo > "${SRA_DIR}/runinfo.csv" || true
  echo "    Wrote ${SRA_DIR}/runinfo.csv (Entrez Direct runinfo) for manual verification."
else
  echo "    Entrez Direct (esearch/efetch) not found -- install with:"
  echo "    conda install -c bioconda entrez-direct"
  echo "    or verify runs manually at:"
  echo "    https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${BIOPROJECT}"
fi

mapfile -t ACCESSIONS < <(yq -r '.sra_accessions[]' "$CONFIG")

for SRR in "${ACCESSIONS[@]}"; do
  echo ">>> prefetch $SRR"
  prefetch --max-size 100G -O "$SRA_DIR" "$SRR"

  echo ">>> fasterq-dump $SRR (split into R1/R2, gzip after)"
  fasterq-dump \
    --split-files \
    --include-technical \
    -e 8 \
    -O "$FASTQ_DIR" \
    "${SRA_DIR}/${SRR}/${SRR}.sra"

  # Galaxy required fastqsanger.gz; STAR/fasterq-dump already emit Sanger-encoded
  # FASTQ, we just need to gzip and name consistently.
  for f in "${FASTQ_DIR}/${SRR}"_*.fastq; do
    [ -e "$f" ] || continue
    gzip -f "$f"
  done
  echo ">>> Done: $SRR -> ${FASTQ_DIR}/${SRR}_1.fastq.gz (+_2, [+_3 if technical/barcode read])"
done

echo ">>> All downloads complete. Files are in ${FASTQ_DIR}"
