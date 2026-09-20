"""Consensus figure script — Reha mouse rescue analysis (manuscript-ready set).

Generates Fig. 2b and Fig. 2c from the intermediate CSVs in `outputs/restoration/`
written by `Reha_rescue_analysis_260919.ipynb`. No re-fit of pseudobulk DESeq2
or pathway GSEA — those upstream computations stay in the notebook; this
script consumes their saved outputs.

  MAIN figures
    1.  F_MAIN_vizA1_DESeq2_log_Reds.pdf
        Per-celltype rescue gene count (log10 horizontal bars, Reds gradient
        encoding rescue rate %). Sort: rescue count desc.
    2.  F_MAIN_pathway_restoration_heatmap_DESeq2.pdf
        Hallmark pathway NES heatmap (left = SED vs Ctrl disease signature,
        right = Ex vs SED restoration), sign-flip & |NES|>1 & FDR<0.10
        filtered, average-linkage hierarchical clustering on disease NES
        rows. Black borders mark rescue (sign-flip + magnitude). `*` = FDR<0.10.

Inputs (all in outputs/restoration/):
  - rescue_summary_per_celltype_DESeq2.csv    (Viz A1)
  - pathway_NES_per_celltype_DESeq2.csv       (pathway heatmap)

Method context kept identical:
  - Pseudobulk DESeq2 (Heumos 2023 standard for between-condition snRNA DE).
  - MI DEG threshold: FDR<0.10, |log2FC|>0.25; rescue classes follow Ma et al. 2020 (Cell).
  - Pathway: Hallmark MSigDB 2020 (locally pinned GMT) via gseapy.prerank,
    100 perm, mouse→human ortholog via mygene.
  - n=2/group: Sham (Ctrl), MI+SED, MI+Ex(Reha6W).
  - Cell types with zero MI DEGs (B cells, Lymphoid Cells) excluded.

Run: `conda activate scanpy && python make_reha_rescue_figures.py`
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

ROOT = Path("/Users/takahiro/Desktop/project/Reha/snRNAseq_scanpy")
OUT = ROOT / "outputs/restoration"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Helvetica", "font.size": 9,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
})

THR_LABEL = "FDR<0.10, |log2FC|>0.25"
QUAD_COLORS = {
    "overshoot": "#7a3f80", "rescue": "#2ca25f",
    "partial": "#9ed99c", "side_effect": "#d9534f", "no_signal": "#cccccc",
}


# ============================================================================
# 1. MAIN — Viz A1 (rescue gene count log bar, Reds rate%)
# ============================================================================
def make_vizA1():
    print("\n[1/2] MAIN — Viz A1 (rescue count log10 bar, Reds rate%)")
    m = pd.read_csv(OUT / "rescue_summary_per_celltype_DESeq2.csv")
    m = m[m["n_MI_DEGs"] > 0].copy()
    m["rescue_pct"] = m["n_rescue"] / m["n_MI_DEGs"] * 100.0
    m = m.set_index("celltype")
    order = m.sort_values("n_rescue", ascending=False).index.tolist()
    m = m.loc[order]

    PCT_MAX = float(np.ceil(max(m["rescue_pct"].max(), 5) / 2.5) * 2.5)
    norm_pct = mpl.colors.Normalize(vmin=0, vmax=PCT_MAX)
    cmap_pct = plt.cm.Reds

    counts = m["n_rescue"].fillna(0).values.astype(float)
    pcts = m["rescue_pct"].fillna(0).values
    n_mis = m["n_MI_DEGs"].fillna(0).astype(int).values
    colors = [cmap_pct(norm_pct(p)) for p in pcts]
    x_log = np.log10(counts + 1)
    max_x = max(x_log.max(), 1)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.5 * len(m))))
    ax.barh(range(len(m)), np.maximum(x_log, max_x * 0.005),
            color=colors, edgecolor="black", linewidth=0.6)
    for i, ct in enumerate(m.index):
        ax.text(max_x * 1.02, i,
                f"{int(counts[i]):>4d}  |  {pcts[i]:>5.1f}%  (of {int(n_mis[i])} MI DEGs)",
                va="center", fontsize=9, family="monospace")
    ax.set_xticks([np.log10(t + 1) for t in [1, 10, 100, 1000]])
    ax.set_xticklabels(["1", "10", "100", "1000"])
    ax.set_xlim(0, max_x * 1.6)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(m.index, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("# rescue genes (log10 scale)")
    ax.set_title(f"MAIN · pseudobulk DESeq2 (Heumos 2023 standard, {THR_LABEL})",
                 loc="left", fontsize=11)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    sm = mpl.cm.ScalarMappable(cmap=cmap_pct, norm=norm_pct)
    fig.colorbar(sm, ax=ax, location="right", shrink=0.7, label="rescue rate %")
    plt.tight_layout()
    out_pdf = FIG / "F_MAIN_vizA1_DESeq2_log_Reds.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"  wrote {out_pdf}")
    return order   # used as celltype column order in the heatmap


# ============================================================================
# 2. MAIN — Pathway restoration heatmap (Hallmark NES, sign-flip + sig)
# ============================================================================
def make_pathway_heatmap(ct_order_from_vizA1):
    print("\n[2/2] MAIN — Pathway restoration heatmap (Hallmark NES)")
    pw = pd.read_csv(OUT / "pathway_NES_per_celltype_DESeq2.csv")
    for c in ["NES", "NOM p-val", "FDR q-val"]:
        if c in pw.columns: pw[c] = pd.to_numeric(pw[c], errors="coerce")

    nes_wide = pw.pivot_table(index=["celltype", "Term"], columns="contrast", values="NES").dropna()
    nes_wide["flip"] = (
        (np.sign(nes_wide["SED_vs_Ctrl"]) != np.sign(nes_wide["Ex_vs_SED"]))
        & (nes_wide["SED_vs_Ctrl"].abs() > 1) & (nes_wide["Ex_vs_SED"].abs() > 1)
    )
    fdr_d = pw[pw.contrast == "SED_vs_Ctrl"].pivot_table(index="Term", columns="celltype", values="FDR q-val")
    fdr_r = pw[pw.contrast == "Ex_vs_SED"].pivot_table(index="Term", columns="celltype", values="FDR q-val")
    sig_either = ((fdr_d < 0.10).any(axis=1)) | ((fdr_r < 0.10).any(axis=1))
    flip_terms = set(nes_wide.reset_index().query("flip").Term.unique().tolist())
    sig_terms = set(sig_either[sig_either].index.tolist())
    pathway_rows = sorted(flip_terms & sig_terms)

    ct_order = [c for c in ct_order_from_vizA1 if c in pw["celltype"].unique()]
    print(f"  pathways with sign-flip & sig (FDR<0.10 in either): {len(pathway_rows)}")
    print(f"  celltypes: {len(ct_order)}")

    if len(pathway_rows) == 0 or len(ct_order) == 0:
        print("  no rows / cols to plot, skipping")
        return

    nes_d = pw[pw.contrast == "SED_vs_Ctrl"].pivot_table(index="Term", columns="celltype", values="NES").reindex(index=pathway_rows, columns=ct_order)
    nes_r = pw[pw.contrast == "Ex_vs_SED"].pivot_table(index="Term", columns="celltype", values="NES").reindex(index=pathway_rows, columns=ct_order)
    pd_d = fdr_d.reindex(index=pathway_rows, columns=ct_order)
    pd_r = fdr_r.reindex(index=pathway_rows, columns=ct_order)

    if len(pathway_rows) > 1:
        Z = linkage(nes_d.fillna(0).values, method="average", metric="euclidean")
        oi = leaves_list(Z)
        pathway_rows = [pathway_rows[i] for i in oi]
        nes_d = nes_d.loc[pathway_rows]; nes_r = nes_r.loc[pathway_rows]
        pd_d = pd_d.loc[pathway_rows]; pd_r = pd_r.loc[pathway_rows]

    abs_max = max(np.nanmax(np.abs(nes_d.values)), np.nanmax(np.abs(nes_r.values)), 2.0)
    norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)

    fig = plt.figure(figsize=(2.5 + 0.55 * len(ct_order) * 2, max(6, 0.32 * len(pathway_rows))))
    gs = GridSpec(1, 3, width_ratios=[len(ct_order), len(ct_order), 0.8], wspace=0.4)
    ax_d = fig.add_subplot(gs[0, 0]); ax_r = fig.add_subplot(gs[0, 1]); ax_cb = fig.add_subplot(gs[0, 2])

    sns.heatmap(nes_d, ax=ax_d, cmap="RdBu_r", norm=norm, cbar=False,
                linewidths=0.4, linecolor="white", mask=nes_d.isna())
    for i, term in enumerate(pathway_rows):
        for j, ct in enumerate(ct_order):
            p = pd_d.iloc[i, j]
            if pd.notna(p) and p < 0.10:
                ax_d.text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=10, color="black")
    ax_d.set_title("SED vs Ctrl  (disease)", fontsize=12)
    ax_d.set_xlabel(""); ax_d.set_ylabel("Hallmark pathway", fontsize=10)
    ax_d.set_xticklabels(ax_d.get_xticklabels(), rotation=45, ha="right")

    sns.heatmap(nes_r, ax=ax_r, cmap="RdBu_r", norm=norm, cbar=False,
                linewidths=0.4, linecolor="white", mask=nes_r.isna())
    for i, term in enumerate(pathway_rows):
        for j, ct in enumerate(ct_order):
            p = pd_r.iloc[i, j]
            if pd.notna(p) and p < 0.10:
                ax_r.text(j + 0.5, i + 0.5, "*", ha="center", va="center", fontsize=10, color="black")
            d_ = nes_d.iloc[i, j]; r_ = nes_r.iloc[i, j]
            if pd.notna(d_) and pd.notna(r_) and (np.sign(d_) != np.sign(r_)) and abs(d_) > 1 and abs(r_) > 1:
                ax_r.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor="black", lw=1.5))
    ax_r.set_title("Ex vs SED  (restoration)", fontsize=12)
    ax_r.set_xlabel(""); ax_r.set_ylabel(""); ax_r.set_yticklabels([])
    ax_r.set_xticklabels(ax_r.get_xticklabels(), rotation=45, ha="right")
    mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.RdBu_r, norm=norm,
                              orientation="vertical", label="NES")
    plt.suptitle(
        "MAIN · Pathway restoration (pseudobulk DESeq2 LFC ranks, Hallmark GSEA, local GMT + mygene ortholog)\n"
        "* = FDR q<0.10 ; black border = sign-flip & |NES|>1 (rescue)",
        fontsize=11, y=1.005,
    )
    out_pdf = FIG / "F_MAIN_pathway_restoration_heatmap_DESeq2.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"  wrote {out_pdf}")


# ============================================================================
def main():
    print(f"Reha rescue figure regeneration → {FIG}")
    ct_order = make_vizA1()
    make_pathway_heatmap(ct_order)
    print("\nAll figures regenerated.")


if __name__ == "__main__":
    main()
