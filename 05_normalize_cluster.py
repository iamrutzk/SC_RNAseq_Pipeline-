#!/usr/bin/env python
"""
Phase 5 (normalize/log/HVG/scale) + Phase 6 (PCA/neighbors/Leiden/UMAP).

Usage:
    python 05_normalize_cluster.py --config config/config.yaml --sample <sample_name>

Input:  results/<sample>/Cells_Filtered_Object.h5ad
Output: results/<sample>/UMAP_Object.h5ad
        results/<sample>/figures/*.png
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
    args = ap.parse_args()

    cfg = load_config(args.config)
    norm_cfg = cfg["normalize"]
    clu_cfg = cfg["cluster"]
    sample = args.sample
    results_dir = os.path.join(cfg["paths"]["results_dir"], sample)
    fig_dir = os.path.join(results_dir, "figures")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(os.path.join(results_dir, "Cells_Filtered_Object.h5ad"))
    print(f">>> Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    # ---- Phase 5.1-5.2: normalize + log1p ----------------------------------
    sc.pp.normalize_total(adata, target_sum=norm_cfg["target_sum"])
    sc.pp.log1p(adata)

    # ---- Phase 5.3: freeze raw (pre-scaling) normalized data --------------
    adata.raw = adata

    # ---- Phase 5.4: HVGs (flag only, keep all genes) -----------------------
    sc.pp.highly_variable_genes(adata, flavor=norm_cfg["hvg_flavor"])
    print(f">>> {adata.var['highly_variable'].sum()} HVGs flagged")

    # ---- Phase 5.5: scale ----------------------------------------------------
    adata_scaled = adata.copy()
    sc.pp.scale(adata_scaled, max_value=norm_cfg["scale_max_value"])

    # ---- Phase 6.1: PCA --------------------------------------------------
    sc.tl.pca(adata_scaled, n_comps=clu_cfg["n_pcs_compute"], svd_solver="arpack")
    sc.pl.pca_variance_ratio(
        adata_scaled, n_pcs=clu_cfg["n_pcs_compute"],
        save=f"_{sample}_elbow.png", show=False,
    )
    # carry PCA coords back onto the log-normalized (non-scaled) object,
    # matching the Galaxy protocol's "use all genes for graph stability"
    adata.obsm["X_pca"] = adata_scaled.obsm["X_pca"]

    # ---- Phase 6.2: neighbors ----------------------------------------------
    sc.pp.neighbors(
        adata,
        n_neighbors=clu_cfg["n_neighbors"],
        n_pcs=clu_cfg["n_pcs_use"],
        metric="euclidean",
    )

    # ---- Phase 6.3: Leiden clustering ---------------------------------------
    sc.tl.leiden(
        adata,
        resolution=clu_cfg["leiden_resolution"],
        random_state=clu_cfg["random_state"],
    )
    print(">>> Cluster sizes:")
    print(adata.obs["leiden"].value_counts())

    # ---- Phase 6.4: UMAP -----------------------------------------------------
    sc.tl.umap(adata, min_dist=clu_cfg["umap_min_dist"])

    # ---- Phase 6.5: plots -----------------------------------------------------
    sc.pl.umap(
        adata, color="leiden", legend_loc="on data", size=50,
        save=f"_{sample}_clusters.png", show=False,
    )

    species = cfg["species"]
    marker_key = f"neutrophil_markers_{species}"
    markers = [g for g in cfg[marker_key] if g in adata.var_names]
    if markers:
        sc.pl.umap(
            adata, color=markers, cmap="viridis",
            save=f"_{sample}_neutrophil_markers.png", show=False,
        )
    else:
        print(f"!! none of the configured {marker_key} genes were found in var_names")

    out_path = os.path.join(results_dir, "UMAP_Object.h5ad")
    adata.write_h5ad(out_path)
    print(f">>> Wrote {out_path}")


if __name__ == "__main__":
    main()
