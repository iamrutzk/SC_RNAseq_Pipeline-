#!/usr/bin/env python
"""
Phase 3.2 (empty-droplet removal) + Phase 4 (QC & filtering), terminal version.

Replaces the Galaxy "DropletUtils" + "Scanpy Inspect and manipulate" +
"Scanpy plot" + "Scanpy filter" tool chain with one script.

Usage:
    python 04_qc_filter.py --config config/config.yaml --sample <sample_name>

Input:
    results/<sample>/starsolo/Solo.out/Gene/raw/{matrix.mtx,barcodes.tsv,features.tsv}
Output:
    results/<sample>/Cells_Filtered_Object.h5ad
    results/<sample>/figures/*.png

--- Optional: statistical empty-droplet removal via DropletUtils::emptyDrops ---
If you specifically want DropletUtils' emptyDrops() test (rather than the
hard thresholds below, which is what your original doc used as its own
documented fallback), call it via rpy2 before this script, e.g.:

    import rpy2.robjects as ro
    ro.r('''
        library(DropletUtils)
        sce <- read10xCounts("results/<sample>/starsolo/Solo.out/Gene/raw")
        out <- emptyDrops(counts(sce))
        keep <- out$FDR <= 0.01 & !is.na(out$FDR)
        sce <- sce[, keep]
        write10xCounts("results/<sample>/dropletutils_filtered", counts(sce))
    ''')

then point this script's --matrix-dir at that filtered folder instead.
"""
import argparse
import os
import yaml
import scanpy as sc

sc.settings.verbosity = 1


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument(
        "--matrix-dir",
        default=None,
        help="override path to raw matrix dir (defaults to STARsolo raw output)",
    )
    args = ap.parse_args()

    cfg = load_config(args.config)
    qc = cfg["qc"]
    sample = args.sample
    results_dir = os.path.join(cfg["paths"]["results_dir"], sample)
    fig_dir = os.path.join(results_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    matrix_dir = args.matrix_dir or os.path.join(
        results_dir, "starsolo", "Solo.out", "Gene", "raw"
    )

    print(f">>> Loading raw matrix from {matrix_dir}")
    adata = sc.read_mtx(os.path.join(matrix_dir, "matrix.mtx")).T
    barcodes = [l.strip() for l in open(os.path.join(matrix_dir, "barcodes.tsv"))]
    features_path = os.path.join(matrix_dir, "features.tsv")
    if not os.path.exists(features_path):
        features_path = os.path.join(matrix_dir, "genes.tsv")
    features = [l.strip().split("\t") for l in open(features_path)]
    adata.obs_names = barcodes
    adata.var_names = [f[1] if len(f) > 1 else f[0] for f in features]
    adata.var["gene_ids"] = [f[0] for f in features]
    adata.var_names_make_unique()
    print(f">>> Loaded: {adata.n_obs} droplets x {adata.n_vars} genes")

    # ---- Phase 3.2 substitute: crude empty-droplet knee filter ------------
    # Keep only barcodes with at least a handful of counts before computing
    # QC metrics -- avoids wasting time/memory on millions of true-empty
    # droplets. This is NOT a replacement for emptyDrops' statistical test,
    # just a practical pre-filter (mirrors DropletUtils' "expected cells"
    # ballpark from the original doc).
    sc.pp.filter_cells(adata, min_counts=1)

    # ---- Phase 4.1: QC metrics --------------------------------------------
    mito_prefix = qc["mito_prefix"]
    adata.var["mito"] = adata.var_names.str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mito"], percent_top=None, log1p=True, inplace=True
    )
    adata.obs.rename(columns={"pct_counts_mito": "pct_counts_mito"}, inplace=True)

    # ---- Phase 4.2: QC visualizations -------------------------------------
    sc.pl.violin(
        adata,
        ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mito"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{sample}_qc_violin.png",
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="log1p_total_counts",
        y="pct_counts_mito",
        save=f"_{sample}_UMIxMito.png",
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="log1p_n_genes_by_counts",
        y="pct_counts_mito",
        save=f"_{sample}_GenesxMito.png",
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts",
        color="pct_counts_mito",
        save=f"_{sample}_GenesxUMI.png",
        show=False,
    )

    n_before = adata.n_obs
    print(f">>> Cells before filtering: {n_before}")

    # ---- Phase 4.3: sequential filtering (same thresholds as your doc) ----
    adata = adata[adata.obs["log1p_n_genes_by_counts"] > qc["min_log1p_n_genes"]]
    adata = adata[adata.obs["log1p_n_genes_by_counts"] < qc["max_log1p_n_genes"]]
    adata = adata[adata.obs["log1p_total_counts"] > qc["min_log1p_total_counts"]]
    adata = adata[adata.obs["log1p_total_counts"] < qc["max_log1p_total_counts"]]
    adata = adata[adata.obs["pct_counts_mito"] < qc["max_pct_mito"]]
    sc.pp.filter_genes(adata, min_cells=qc["min_cells_per_gene"])

    print(f">>> Cells after filtering: {adata.n_obs} ({adata.n_obs/n_before:.1%} kept)")
    print(f">>> Genes after filtering: {adata.n_vars}")

    out_path = os.path.join(results_dir, "Cells_Filtered_Object.h5ad")
    adata.write_h5ad(out_path)
    print(f">>> Wrote {out_path}")


if __name__ == "__main__":
    main()
