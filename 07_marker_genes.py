#!/usr/bin/env python
"""
Phase 8: cluster marker genes (rank_genes_groups), heatmap, exported table.

Usage:
    python 07_marker_genes.py --config config/config.yaml --sample <sample_name> \
        [--input Neutrophil_Subset.h5ad]   # defaults to UMAP_Object.h5ad (all cells)

Output: results/<sample>/markers_top<N>_per_cluster.csv
        results/<sample>/figures/*_heatmap.png
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
    ap.add_argument("--input", default="UMAP_Object.h5ad")
    args = ap.parse_args()

    cfg = load_config(args.config)
    mk_cfg = cfg["markers"]
    sample = args.sample
    results_dir = os.path.join(cfg["paths"]["results_dir"], sample)
    fig_dir = os.path.join(results_dir, "figures")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(os.path.join(results_dir, args.input))

    # ---- Phase 8.1: rank_genes_groups ---------------------------------------
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method=mk_cfg["method"],
        n_genes=mk_cfg["n_genes"],
        pts=True,
    )

    # ---- Phase 8.2: heatmap of top 10 per cluster ---------------------------
    top10 = {}
    for cluster in adata.obs["leiden"].cat.categories:
        names = adata.uns["rank_genes_groups"]["names"][cluster][:10]
        top10[cluster] = list(names)
    all_top10_genes = sorted(set(g for genes in top10.values() for g in genes))

    sc.pl.heatmap(
        adata,
        var_names=all_top10_genes,
        groupby="leiden",
        cmap="Reds",
        save=f"_{sample}_marker_heatmap.png",
        show=False,
    )

    # ---- Phase 8.3: export full ranked table --------------------------------
    records = []
    for cluster in adata.obs["leiden"].cat.categories:
        names = adata.uns["rank_genes_groups"]["names"][cluster]
        lfc = adata.uns["rank_genes_groups"]["logfoldchanges"][cluster]
        pvals_adj = adata.uns["rank_genes_groups"]["pvals_adj"][cluster]
        for gene, l, p in zip(names, lfc, pvals_adj):
            records.append(
                {"cluster": cluster, "gene": gene, "log2FC": l, "padj": p}
            )
    df = pd.DataFrame(records)
    out_csv = os.path.join(
        results_dir, f"markers_top{mk_cfg['n_genes']}_per_cluster.csv"
    )
    df.to_csv(out_csv, index=False)
    print(f">>> Wrote {out_csv} ({len(df)} rows)")


if __name__ == "__main__":
    main()
