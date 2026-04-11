#!/usr/bin/env python3
# ============================================================
# Fig18BD_Velocity.py
# Figure 18B/D: RNA Velocity (scVelo Stochastic, Lancet palette)
# Input:  ./out/velocyto/para_velocity.h5ad
# Output: ./figures/Fig18B_*.pdf, Fig18D_*.pdf, Fig18BD_*.pdf
#         ./out/velocyto/velocity_metadata.csv
# Usage:  conda activate scanpy-env && python scripts/velocyto/Fig18BD_Velocity.py
# ============================================================

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.backends.backend_pdf as _pdf
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import scvelo as scv

scv.settings.verbosity = 3
# --- Monkey-patch: PDF backend non-finite safety ---
_orig_pdfRepr = _pdf.pdfRepr
def _safe_pdfRepr(obj):
    if isinstance(obj, (float, np.floating)):
        if not np.isfinite(obj):
            obj = 0.0
    elif isinstance(obj, np.ndarray):
        obj = obj.copy()
        obj[~np.isfinite(obj)] = 0.0
    return _orig_pdfRepr(obj)
_pdf.pdfRepr = _safe_pdfRepr

# --- Paths ---
basedir = "./out/velocyto/"
figdir  = "./figures/"
os.makedirs(figdir, exist_ok=True)

# --- Lancet palette (ggsci, 8 colors) ---
LANCET_8 = [
    "#00468B", "#ED0000", "#42B540", "#0099B4",
    "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6",
]

def assign_lancet_colors(adata, col):
    cats = list(adata.obs[col].cat.categories)
    adata.uns[f"{col}_colors"] = [LANCET_8[i % 8] for i in range(len(cats))]
# ============================================================
# Step 1: Load h5ad & preprocess
# ============================================================
adata = scv.read(os.path.join(basedir, "para_velocity.h5ad"))
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

for col in ["celltype", "group", "orig.ident", "seurat_clusters"]:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

for col in adata.obs.select_dtypes(include=[np.floating]).columns:
    mask = ~np.isfinite(adata.obs[col].values)
    if mask.any():
        adata.obs.loc[mask, col] = 0.0

assign_lancet_colors(adata, "celltype")
assign_lancet_colors(adata, "group")
# ============================================================
# Step 2: scVelo stochastic velocity
# ============================================================
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
if "X_pca" in adata.obsm:
    scv.pp.moments(adata, n_pcs=30, use_rep="X_pca", n_neighbors=30)
else:
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)
print(f"Velocity genes: {adata.var['velocity_genes'].sum()} / {adata.n_vars}")
# --- Global UMAP range (per-group subplots use same axes) ---
umap_coords = adata.obsm["X_umap"]
x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
pad_x = (x_max - x_min) * 0.05
pad_y = (y_max - y_min) * 0.05
groups = sorted(adata.obs["group"].unique().tolist())
n_groups = len(groups)

# ============================================================
# Figure 18B: Velocity stream (celltype & group)
# ============================================================
for color_key in ["celltype", "group"]:
    fig, ax = plt.subplots(figsize=(8, 7))
    scv.pl.velocity_embedding_stream(
        adata, basis="umap", color=color_key,
        title=f"RNA Velocity — {color_key}",
        ax=ax, show=False, legend_loc="right margin")
    plt.tight_layout()
    fig.savefig(os.path.join(figdir, f"Fig18B_stream_{color_key}.pdf"),
                bbox_inches="tight")
    plt.close("all")
# --- Combined 3-panel ---
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="celltype",
    title="celltype", ax=axes[0], show=False, legend_loc="right margin")
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="group",
    title="group", ax=axes[1], show=False, legend_loc="right margin")
sc.pl.umap(adata, color="velocity_pseudotime",
    ax=axes[2], show=False, title="pseudotime")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18B_stream_combined.pdf"),
            bbox_inches="tight")
plt.close("all")
# ============================================================
# Figure 18D: Per-group velocity stream
# ============================================================
fig, axes = plt.subplots(1, n_groups, figsize=(8 * n_groups, 7))
if n_groups == 1:
    axes = [axes]

# 全局 celltype → color 映射字典（按 categories 顺序）
_global_ct_cats   = list(adata.obs["celltype"].cat.categories)
_global_ct_colors = list(adata.uns["celltype_colors"])
_ct_color_map     = dict(zip(_global_ct_cats, _global_ct_colors))

for i, g in enumerate(groups):
    adata_sub = adata[adata.obs["group"] == g].copy()
    # 移除子集中不存在的 categories，再按全局映射赋色
    adata_sub.obs["celltype"] = adata_sub.obs["celltype"].cat.remove_unused_categories()
    sub_cats = list(adata_sub.obs["celltype"].cat.categories)
    adata_sub.uns["celltype_colors"] = [_ct_color_map[c] for c in sub_cats]
    sc.pp.neighbors(adata_sub, n_neighbors=30, use_rep="X_pca")
    scv.pp.moments(adata_sub)
    scv.tl.velocity_graph(adata_sub)
    scv.pl.velocity_embedding_stream(
        adata_sub, basis="umap", color="celltype",
        title=g, ax=axes[i], show=False,
        legend_loc="on data", legend_fontsize=16)
    for txt in axes[i].texts:
        txt.set_color("white")
        txt.set_fontweight("bold")
        txt.set_fontsize(16)
        txt.set_path_effects([
            pe.withStroke(linewidth=2.5, foreground="black"),
            pe.Normal()])
    axes[i].set_xlim(x_min - pad_x, x_max + pad_x)
    axes[i].set_ylim(y_min - pad_y, y_max + pad_y)
    axes[i].set_aspect("equal", adjustable="box")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18D_stream_pergroup.pdf"),
            bbox_inches="tight")
plt.close("all")

# ============================================================
# Arrow embedding (global + per-group)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
scv.pl.velocity_embedding(
    adata, basis="umap", color="celltype",
    arrow_length=3, arrow_size=2, title="Arrows — celltype",
    ax=axes[0], show=False, legend_loc="right margin")
scv.pl.velocity_embedding(
    adata, basis="umap", color="group",
    arrow_length=3, arrow_size=2, title="Arrows — group",
    ax=axes[1], show=False, legend_loc="right margin")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18BD_arrows.pdf"),
            bbox_inches="tight")
plt.close("all")
fig, axes = plt.subplots(1, n_groups, figsize=(8 * n_groups, 7))
if n_groups == 1:
    axes = [axes]
for i, g in enumerate(groups):
    adata_sub = adata[adata.obs["group"] == g].copy()
    # 移除子集中不存在的 categories，再按全局映射赋色
    adata_sub.obs["celltype"] = adata_sub.obs["celltype"].cat.remove_unused_categories()
    sub_cats = list(adata_sub.obs["celltype"].cat.categories)
    adata_sub.uns["celltype_colors"] = [_ct_color_map[c] for c in sub_cats]
    sc.pp.neighbors(adata_sub, n_neighbors=30, use_rep="X_pca")
    scv.pp.moments(adata_sub)
    scv.tl.velocity_graph(adata_sub)
    scv.pl.velocity_embedding(
        adata_sub, basis="umap", color="celltype",
        arrow_length=3, arrow_size=2, title=g,
        ax=axes[i], show=False, legend_loc="on data", legend_fontsize=16)
    # 白色字体 + 黑色描边（与 stream 子图一致）
    for txt in axes[i].texts:
        txt.set_color("white")
        txt.set_fontweight("bold")
        txt.set_fontsize(16)
        txt.set_path_effects([
            pe.withStroke(linewidth=2.5, foreground="black"),
            pe.Normal()])
    axes[i].set_xlim(x_min - pad_x, x_max + pad_x)
    axes[i].set_ylim(y_min - pad_y, y_max + pad_y)
    axes[i].set_aspect("equal", adjustable="box")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18BD_arrows_pergroup.pdf"),
            bbox_inches="tight")
plt.close("all")
# ============================================================
# Confidence & pseudotime
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="velocity_confidence",
    ax=axes[0], show=False, title="Velocity Confidence")
sc.pl.umap(adata, color="velocity_length",
    ax=axes[1], show=False, title="Velocity Length")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18BD_confidence_length.pdf"),
            bbox_inches="tight")
plt.close("all")

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color="velocity_pseudotime",
    ax=axes[0], show=False, title="Velocity Pseudotime", color_map="gnuplot")
sc.pl.umap(adata, color="celltype",
    ax=axes[1], show=False, title="celltype reference")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18BD_pseudotime.pdf"),
            bbox_inches="tight")
plt.close("all")
fig, ax = plt.subplots(figsize=(8, 5))
sc.pl.violin(adata, keys="velocity_pseudotime", groupby="celltype",
    rotation=45, ax=ax, show=False)
ax.set_title("Velocity Pseudotime by celltype")
plt.tight_layout()
fig.savefig(os.path.join(figdir, "Fig18BD_pseudotime_violin.pdf"),
            bbox_inches="tight")
plt.close("all")

# ============================================================
# Rank velocity genes & save
# ============================================================
scv.tl.rank_velocity_genes(adata, groupby="celltype", min_corr=0.3)
df_genes = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]).head(10)
df_genes.to_csv(os.path.join(basedir, "rank_velocity_genes.csv"), index=False)

meta = adata.obs[
    ["celltype", "group", "velocity_pseudotime",
     "velocity_confidence", "velocity_length"]].copy()
meta.to_csv(os.path.join(basedir, "velocity_metadata.csv"))
adata.write(os.path.join(basedir, "para_velocity_result.h5ad"))
print("=== DONE ===")