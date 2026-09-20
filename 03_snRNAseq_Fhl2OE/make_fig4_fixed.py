#!/usr/bin/env python
"""Assemble the Fig. 4e–h panels from the corrected Fhl2-OE run.

Same panel roles and visual style as the published figure:
  e  UMAP of cardiomyocytes coloured by condition (SED / Ex / Fhl2 OE)
  f  gene-score violins (cardiac contraction, metabolites-energy, inflammation)
  g  single-gene violins (Nppb, Myh7, Tlr4)
  h  GO/pathway enrichment for genes up in Fhl2 OE vs Ex

Significance is computed from the data (Kruskal-Wallis across the three groups, then
pairwise two-sided Mann-Whitney U on the annotated pairs) — never assumed.

Run with the `scanpy` env:
  ~/miniforge3/envs/scanpy/bin/python make_fig4_fixed.py
"""
import warnings
import numpy as np, pandas as pd, scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import kruskal, mannwhitneyu

warnings.filterwarnings("ignore")
BASE = "/Users/takahiro/Desktop/project/Reha/snRNAseq_scanpy/Fhl2OE"
AD   = f"{BASE}/adata_fixed/adata_Fhl2_analysed_fixed.h5ad"
FIG  = f"{BASE}/figures_fixed"

FHL2 = r"$\it{Fhl2}$ OE"
GROUPS = ["SED", "Ex", FHL2]
LABELS = ["SED", "Ex", "Fhl2 OE"]
COLORS = ["#2ca02c", "#ff7f0e", "#d62728"]      # as in the original notebook
PAIRS  = [("SED", "Ex"), ("SED", FHL2)]
GENE_PANEL = ["Nppb", "Cox4i1", "Tlr4"]   # Fig. 4g, as listed in the legend of the 10 Sep 26 draft

mpl.rcParams.update({"font.size": 9, "axes.grid": False, "pdf.fonttype": 42})


def stars(p):
    return "****" if p <= 1e-4 else "***" if p <= 1e-3 else "**" if p <= 1e-2 else "*" if p <= 5e-2 else "ns"


def violin_panel(ax, adata, feature, from_obs, title, ylabel=None, italic_title=False):
    if from_obs:
        vals = adata.obs[feature].astype(float).values
    else:
        vals = np.asarray(adata.raw[:, feature].X.todense()).ravel()
    grp = adata.obs["type"].astype(str).values
    data = [vals[grp == g] for g in GROUPS]

    H, p_kw = kruskal(*data)
    parts = ax.violinplot(data, positions=range(3), showextrema=False, widths=0.85)
    for pc, c in zip(parts["bodies"], COLORS):
        pc.set_facecolor(c); pc.set_alpha(0.85); pc.set_edgecolor("black"); pc.set_linewidth(0.4)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        idx = rng.choice(len(d), size=min(len(d), 400), replace=False)
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(idx)), d[idx], s=1.2,
                   color="black", alpha=0.25, linewidths=0, zorder=3)

    lo = min(d.min() for d in data); hi = max(d.max() for d in data)
    span = hi - lo
    y = hi + span * 0.06
    res = {"feature": feature, "kruskal_H": H, "kruskal_p": p_kw}
    for k, (a, b) in enumerate(PAIRS):
        ia, ib = GROUPS.index(a), GROUPS.index(b)
        U, p = mannwhitneyu(data[ia], data[ib], alternative="two-sided")
        res[f"p_{LABELS[ia]}_vs_{LABELS[ib]}"] = p
        ax.plot([ia, ia, ib, ib], [y, y + span*0.03, y + span*0.03, y], lw=0.8, c="black")
        ax.text((ia + ib) / 2, y + span*0.035, stars(p), ha="center", va="bottom", fontsize=10)
        y += span * 0.14
    ax.set_xticks(range(3)); ax.set_xticklabels(["SED", "Ex", ""], fontsize=8)
    tr = ax.get_xaxis_transform()
    ax.text(2 - 0.06, -0.035, "Fhl2", ha="right", va="top", style="italic",
            fontsize=8, transform=tr)
    ax.text(2 - 0.06, -0.035, " OE", ha="left", va="top", fontsize=8, transform=tr)
    ax.set_ylim(lo - span * 0.05, y + span * 0.05)
    ax.set_ylabel(ylabel or feature.replace("_", " "), fontsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold",
                 style="italic" if italic_title else "normal")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    return res


def main():
    adata = sc.read_h5ad(AD)
    n = adata.obs["type"].value_counts()
    print("cells:", dict(n))

    fig = plt.figure(figsize=(15, 8), dpi=150)
    gs = GridSpec(2, 4, figure=fig, hspace=0.42, wspace=0.38)

    # ---- e : UMAP -------------------------------------------------------
    axe = fig.add_subplot(gs[0, 0])
    emb = adata.obsm["X_umap"]; grp = adata.obs["type"].astype(str).values
    for g, c, lab in zip(GROUPS, COLORS, LABELS):
        m = grp == g
        axe.scatter(emb[m, 0], emb[m, 1], s=5, c=c, alpha=0.8, linewidths=0,
                    label=f"{lab} ({m.sum():,})")
    axe.set_xticks([]); axe.set_yticks([])
    axe.set_xlabel("UMAP1", fontsize=8); axe.set_ylabel("UMAP2", fontsize=8)
    axe.legend(markerscale=6, fontsize=7, frameon=False, loc="best")
    axe.set_title("Cardiomyocyte snRNA-seq", fontsize=9, fontweight="bold", loc="left")
    for s in ("top", "right"): axe.spines[s].set_visible(False)

    # ---- f : gene scores ------------------------------------------------
    rows = []
    score_specs = [("cardiac_contraction_score", "Cardiac contraction"),
                   ("metabolites_energy_score",  "Metabolites–energy"),
                   ("inflammation_score",        "Inflammation")]
    for j, (feat, title) in enumerate(score_specs):
        ax = fig.add_subplot(gs[0, j + 1])
        rows.append(violin_panel(ax, adata, feat, True,
                                 title, "gene score"))

    # ---- g : single genes ----------------------------------------------
    for j, gene in enumerate(GENE_PANEL):
        ax = fig.add_subplot(gs[1, j + 4 - len(GENE_PANEL)])
        rows.append(violin_panel(ax, adata, gene, False, gene, "expression",
                                 italic_title=True))

    # ---- h : enrichment -------------------------------------------------
    axh = fig.add_subplot(gs[1, 0:4 - len(GENE_PANEL)])
    try:
        e = pd.read_csv(f"{BASE}/adata_fixed/enrichr_full_up_in_Fhl2OE.csv")
        e = e[e["Gene_set"] == "GO_Biological_Process_2023"]
        e = e[e["Adjusted P-value"] < 0.05].copy()
        e["-log10(adjusted P-value)"] = -np.log10(e["Adjusted P-value"])
        e = e.sort_values("-log10(adjusted P-value)", ascending=False)
        top = e.head(10)
        axh.barh(width=top["-log10(adjusted P-value)"],
                 y=[x.split(" (")[0] for x in top["Term"]], color="darkred")
        axh.invert_yaxis()
        axh.set_xlabel("-log10(adjusted P-value)", fontsize=8)
        axh.tick_params(axis="y", labelsize=6)
        axh.set_title("GO_Biological_Process_2023", fontsize=9, fontweight="bold", loc="left")
        for s in ("top", "right"): axh.spines[s].set_visible(False)
    except FileNotFoundError:
        axh.axis("off")

    fig.suptitle("Figure 4e–h (corrected run, 2026-09-08): "
                 f"SED {n.get('SED',0):,} / Ex {n.get('Ex',0):,} / Fhl2 OE {n.get(FHL2,0):,} nuclei",
                 fontsize=10)
    fig.savefig(f"{FIG}/Fig4_corrected_panels.png", bbox_inches="tight", dpi=200)
    fig.savefig(f"{FIG}/Fig4_corrected_panels.pdf", bbox_inches="tight")
    stats = pd.DataFrame(rows)
    stats.to_csv(f"{BASE}/adata_fixed/Fig4_statistics_fixed.csv", index=False)
    print(stats.to_string(index=False))
    print("saved:", f"{FIG}/Fig4_corrected_panels.png")


if __name__ == "__main__":
    main()
