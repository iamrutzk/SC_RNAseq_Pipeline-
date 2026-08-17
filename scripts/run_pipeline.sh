#!/usr/bin/env bash
# Runs the full pipeline end-to-end for one sample.
# Usage: run_pipeline.sh config/config.yaml <sample_name> \
#          <cDNA_R2_1.fastq.gz[,R2_lane2...]> <barcode_R1_1.fastq.gz[,R1_lane2...]> \
#          [neutrophil_cluster_ids_csv]
#
# Requires `yq` (https://github.com/mikefarah/yq) for parsing config.yaml.
set -euo pipefail

CONFIG="${1:?config path, e.g. config/config.yaml}"
SAMPLE="${2:?sample name}"
CDNA_READS="${3:?comma-separated cDNA fastq.gz files}"
BARCODE_READS="${4:?comma-separated barcode fastq.gz files}"
NEUT_CLUSTERS="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==================================================================="
echo " Sample: $SAMPLE"
echo " Config: $CONFIG"
echo "==================================================================="

echo "### [1/6] Reference genome + STAR index"
bash "${SCRIPT_DIR}/02_get_reference.sh" "$CONFIG"

echo "### [2/6] STARsolo alignment + counting"
bash "${SCRIPT_DIR}/03_run_starsolo.sh" "$CONFIG" "$SAMPLE" "$CDNA_READS" "$BARCODE_READS"

echo "### [3/6] QC metrics + sequential filtering"
python "${SCRIPT_DIR}/04_qc_filter.py" --config "$CONFIG" --sample "$SAMPLE"

echo "### [4/6] Normalize, HVG, PCA, Leiden clustering, UMAP"
python "${SCRIPT_DIR}/05_normalize_cluster.py" --config "$CONFIG" --sample "$SAMPLE"

echo "### [5/6] Neutrophil identification + NETosis scoring"
if [ -n "$NEUT_CLUSTERS" ]; then
  python "${SCRIPT_DIR}/06_netosis_scoring.py" --config "$CONFIG" --sample "$SAMPLE" \
    --neutrophil-clusters "$NEUT_CLUSTERS"
else
  echo ">>> No cluster IDs given yet -- printing marker table only."
  echo ">>> Inspect results/${SAMPLE}/figures/*_neutrophil_markers.png and"
  echo ">>> results/${SAMPLE}/cluster_neutrophil_marker_means.csv, then re-run:"
  echo "    bash ${SCRIPT_DIR}/run_pipeline.sh $CONFIG $SAMPLE '$CDNA_READS' '$BARCODE_READS' <cluster_ids>"
  python "${SCRIPT_DIR}/06_netosis_scoring.py" --config "$CONFIG" --sample "$SAMPLE"
  exit 0
fi

echo "### [6/6] Cluster marker genes (differential expression)"
python "${SCRIPT_DIR}/07_marker_genes.py" --config "$CONFIG" --sample "$SAMPLE" --input UMAP_Object.h5ad
python "${SCRIPT_DIR}/07_marker_genes.py" --config "$CONFIG" --sample "$SAMPLE" --input Neutrophil_Subset.h5ad

echo "==================================================================="
echo " DONE. Key outputs in results/${SAMPLE}/:"
echo "   Neutrophil_NETosis_Scored.h5ad"
echo "   markers_top*_per_cluster.csv"
echo "   figures/*.png"
echo "==================================================================="
