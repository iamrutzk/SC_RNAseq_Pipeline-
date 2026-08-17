#!/usr/bin/env python
"""
Phase 7: identify neutrophil clusters, subset them, score the NETosis
gene signature. Runs natively in Python -- no Jupyter round-trip / no
re-upload to Galaxy needed.

Usage:
    python 06_netosis_scoring.py --config config/config.yaml --sample <sample_name> \
        --neutrophil-clusters 2,5,8

If --neutrophil-clusters is omitted, the script prints per-cluster mean
expression of the neutrophil marker genes so you can pick clusters
yourself, then exits without subsetting (run again with the flag set).

Input:  results/<sample>/UMAP_Object.h5ad
Output: results/<sample>/Neutrophil_Subset.h5ad
        results/<sample>/Neutrophil_NETosis_Scored.h5ad
        results/<sample>/figures/*.png
"""
import argparse
import os
import yaml
import scanpy as sc
import pandas as pd

sc.settings.verbosity = 1


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument(
        "--neutrophil-clusters",
        default=None,
        help="comma-separated leiden cluster IDs identified as neutrophils, e.g. 2,5,8",
    )
    args = ap.parse_args()

    cfg = load_config(args.config)
    species = cfg["species"]
    sample = args.sample
    results_dir = os.path.join(cfg["paths"]["results_dir"], sample)
    fig_dir = os.path.join(results_dir, "figures")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(os.path.join(results_dir, "UMAP_Object.h5ad"))

    marker_key = f"neutrophil_markers_{species}"
    markers = [g for g in cfg[marker_key] if g in adata.var_names]

    # ---- Phase 7.1: help identify neutrophil clusters -----------------------
    mean_expr = (
        sc.get.obs_df(adata, keys=markers + ["leiden"])
        .groupby("leiden")
        .mean(numeric_only=True)
    )
    print(">>> Mean neutrophil-marker expression per cluster:")
    print(mean_expr.to_string())
    mean_expr.to_csv(os.path.join(results_dir, "cluster_neutrophil_marker_means.csv"))

    if not args.neutrophil_clusters:
        print(
            "\n>>> No --neutrophil-clusters given. Inspect the table above "
            "(and figures/*_neutrophil_markers.png) then re-run with e.g.:\n"
            f"    python {__file__} --config {args.config} --sample {sample} "
            "--neutrophil-clusters 2,5,8"
        )
        return

    clusters = [c.strip() for c in args.neutrophil_clusters.split(",")]

    # ---- Phase 7.2: subset neutrophils --------------------------------------
    neut = adata[adata.obs["leiden"].isin(clusters)].copy()
    print(f">>> Subset to {neut.n_obs} neutrophil cells (clusters {clusters})")
    neut.write_h5ad(os.path.join(results_dir, "Neutrophil_Subset.h5ad"))

    # ---- Phase 7.3: NETosis signature scoring -------------------------------
    netosis_key = f"netosis_genes_{species}"
    netosis_genes = [g for g in cfg[netosis_key] if g in neut.var_names]
    missing = set(cfg[netosis_key]) - set(netosis_genes)
    if missing:
        print(f"!! NETosis genes not found in dataset (skipped): {sorted(missing)}")
    if not netosis_genes:
        raise SystemExit(
            "No NETosis signature genes found in var_names -- check gene symbol "
            "casing (mouse=Title case, human=UPPER case) and re-check config.yaml."
        )

    sc.tl.score_genes(neut, gene_list=netosis_genes, score_name="NETosis_Score")

    # ---- Phase 7.4: visualize -----------------------------------------------
    sc.pl.umap(
        neut, color="NETosis_Score", cmap="RdYlBu_r", legend_loc="on data",
        save=f"_{sample}_NETosis_Score.png", show=False,
    )

    out_path = os.path.join(results_dir, "Neutrophil_NETosis_Scored.h5ad")
    neut.write_h5ad(out_path)
    print(f">>> Wrote {out_path}")
    print(
        neut.obs["NETosis_Score"]
        .describe()
        .to_string()
    )


if __name__ == "__main__":
    main()
