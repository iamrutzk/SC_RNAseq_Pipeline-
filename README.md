# Tumor-Infiltrating Neutrophil / NETosis scRNA-seq Pipeline (Command-Line, No Galaxy)

This project reimplements the Galaxy workflow you were given as a pure
Linux command-line pipeline: `sra-tools` → `STAR/STARsolo` → `Python/Scanpy`.
Everything runs from bash + Python scripts, no web UI required.

## 1. Dataset — what's actually on NCBI SRA

I checked the accessions referenced in your original protocol against
NCBI GEO/SRA directly:

| Accession | What it is | Raw FASTQ on SRA? |
|---|---|---|
| **GSE176078** (Wu et al. 2021, breast cancer atlas, 26 tumors) | scSubtype breast cancer atlas | **No.** Raw human reads were deposited in **EGA** (`EGAS00001005173`), a *controlled-access* repository, for patient-privacy reasons. Only processed count matrices are on GEO. You cannot `prefetch` FASTQ for this one. |
| **SRP299847 / PRJEB43259** | Generic placeholders in the original doc | Not a real, verified single dataset for this project — treat these as examples to replace, not literal accessions to use. |
| **GSE127465** (Zilionis et al. 2019, *Immunity*) — "Single-cell transcriptomics of human and mouse lung cancers reveals conserved myeloid populations" | Tumor-infiltrating myeloid cells incl. **neutrophils**, human NSCLC + matched mouse lung tumor model (KP1.9), CD45+ / RBC-depleted sorted cells | **Partially.** Human raw reads are also privacy-restricted and *not* on SRA. **Mouse raw FASTQ *is* public** on SRA under BioProject **`PRJNA524857`** (SRA study `SRP187083`) — 2 healthy mice + 2 tumor-bearing mice, CD45+ sorted lung cells, including a well-characterized neutrophil compartment. |

**Recommendation:** use the **mouse arm of GSE127465** (`PRJNA524857`) as your
real, downloadable, tumor-infiltrating-neutrophil dataset. It is directly
on-topic (myeloid/neutrophil states in a lung tumor model), and it's one of
the few tumor-neutrophil scRNA-seq studies where raw reads are actually
open-access on SRA rather than locked behind EGA/dbGaP.

⚠️ **Important chemistry note:** this dataset was generated with **inDrop**,
not 10x Chromium. inDrop barcode/UMI structure (barcode split length,
linker sequence, UMI length) differs by inDrop version and is **not**
the same as the 10x parameters in your original Galaxy doc. Do not reuse
the 10x `--soloCBstart/len`/`--soloUMIstart/len` defaults blindly — check
the actual read structure first (`scripts/02b_inspect_read_structure.sh`)
and consult the Klein-lab inDrops docs before finalizing STARsolo barcode
parameters. If you'd rather avoid this complication entirely, swap in any
10x Chromium tumor-neutrophil dataset you have access to — the rest of the
pipeline (from FASTQ onward) is chemistry-agnostic except for that one step.

If your actual goal is GSE176078 (breast cancer) specifically, you have two
options instead of SRA: (1) apply for EGA controlled access
(`EGAS00001005173`), or (2) skip realignment entirely and start straight
from GEO's processed matrices
(`GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz`), which drops you in at
Phase 4 (QC) of this pipeline.

## 2. Requirements

```bash
# Conda/mamba environment covering the whole pipeline
mamba env create -f config/environment.yml
conda activate netosis-scrnaseq
```

Installs: `sra-tools`, `star` (STARsolo), `samtools`, `scanpy`,
`anndata`, `leidenalg`, `python-igraph`, `scikit-misc`.

Disk/compute: reference genome + index ~30 GB, raw FASTQ ~5–15 GB per
sample, STAR index build needs ~32 GB RAM. Run on a workstation/cluster
node, not a laptop.

## 3. Pipeline (mirrors the Galaxy phases 1:1)

```
scripts/
  01_download_sra.sh          # Phase 2: prefetch + fasterq-dump FASTQ from SRA
  02_get_reference.sh         # Phase 1/3: download genome FASTA + GTF, build STAR index
  02b_inspect_read_structure.sh  # sanity-check R1/R2 length & barcode structure before aligning
  03_run_starsolo.sh          # Phase 3.1: STARsolo alignment + counting -> raw matrix
  04_qc_filter.py             # Phase 3.2 + 4: QC metrics, plots, sequential filtering
  05_normalize_cluster.py     # Phase 5+6: normalize, log1p, HVG, scale, PCA, neighbors, Leiden, UMAP
  06_netosis_scoring.py       # Phase 7: neutrophil marker check + NETosis module score
  07_marker_genes.py          # Phase 8: rank_genes_groups + heatmap + marker table export
  run_pipeline.sh             # runs 01-07 in order for one sample
config/
  config.yaml                 # every threshold/parameter in one place (edit this, not the scripts)
  environment.yml
```

## 4. Quick start

```bash
cd neutrophil-netosis-scrnaseq
mamba env create -f config/environment.yml && conda activate netosis-scrnaseq

# edit config/config.yaml: set SRA accessions, species (mouse/GRCm39 by default), paths

bash scripts/run_pipeline.sh config/config.yaml
```

Final output: `results/<sample>/Neutrophil_NETosis_Scored.h5ad`,
QC plots in `results/<sample>/figures/`, marker gene table in
`results/<sample>/markers_top50_per_cluster.csv`.

## 5. Where this differs from the Galaxy protocol

- DropletUtils (R/Bioconductor) is not used; empty-droplet removal is done
  with the same hard thresholds your doc already specified as a fallback
  (>300 genes, >500 UMIs, <4.5% mito) applied directly in Scanpy. If you
  want DropletUtils' `emptyDrops` statistical test specifically, there's a
  commented-out R block in `04_qc_filter.py`'s docstring showing how to
  call it via `rpy2` instead.
- NETosis scoring uses `scanpy.tl.score_genes` exactly as your doc's
  Jupyter fallback did — no Galaxy re-upload step needed since everything
  stays in one Python process.
- Cell-cell communication (Phase 9) is left as an optional add-on
  (`squidpy`/`liana-py` are pure-Python CellChat alternatives that don't
  need a notebook round-trip) — ask if you want that script too.
