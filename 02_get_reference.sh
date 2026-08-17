#!/usr/bin/env bash
# Phase 1/3 (terminal version): download genome FASTA + GTF, build STAR index.
# Usage: 02_get_reference.sh config/config.yaml
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
SPECIES=$(yq -r '.species' "$CONFIG")
REF_DIR=$(yq -r '.paths.ref_dir' "$CONFIG")
FASTA_URL=$(yq -r ".reference.${SPECIES}.fasta_url" "$CONFIG")
GTF_URL=$(yq -r ".reference.${SPECIES}.gtf_url" "$CONFIG")
OVERHANG=$(yq -r '.reference.star_index_sjdb_overhang' "$CONFIG")
THREADS=$(yq -r '.reference.star_index_threads' "$CONFIG")

mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo ">>> Downloading genome FASTA ($SPECIES)"
FASTA_GZ=$(basename "$FASTA_URL")
[ -f "${FASTA_GZ%.gz}" ] || { curl -L -o "$FASTA_GZ" "$FASTA_URL" && gunzip -f "$FASTA_GZ"; }

echo ">>> Downloading GTF annotation ($SPECIES)"
GTF_GZ=$(basename "$GTF_URL")
[ -f "${GTF_GZ%.gz}" ] || { curl -L -o "$GTF_GZ" "$GTF_URL" && gunzip -f "$GTF_GZ"; }

FASTA="${FASTA_GZ%.gz}"
GTF="${GTF_GZ%.gz}"

echo ">>> Building STAR index (this needs ~32GB RAM and 30-60 min)"
mkdir -p star_index
STAR \
  --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$OVERHANG" \
  --runThreadN "$THREADS"

echo ">>> STAR index ready at ${REF_DIR}/star_index"

# 10x barcode whitelists (only needed if chemistry is 10x_v2/10x_v3)
CHEMISTRY=$(yq -r '.chemistry' "$CONFIG")
if [[ "$CHEMISTRY" == 10x_* ]]; then
  WL_URL=$(yq -r ".solo_${CHEMISTRY}.whitelist_url" "$CONFIG")
  WL_FILE=$(basename "$WL_URL")
  echo ">>> Downloading 10x barcode whitelist ($CHEMISTRY)"
  curl -L -o "$WL_FILE" "$WL_URL"
  [[ "$WL_FILE" == *.gz ]] && gunzip -f "$WL_FILE"
  echo ">>> Whitelist ready: ${REF_DIR}/${WL_FILE%.gz}"
fi
