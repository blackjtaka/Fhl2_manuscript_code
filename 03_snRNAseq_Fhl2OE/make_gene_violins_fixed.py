#!/usr/bin/env python
"""Violin plots for every gene the original notebook showed, corrected run.

Same style as the Fig. 4 panels: violin + strip, brackets over SED-Ex and SED-Fhl2 OE,
stars from Kruskal-Wallis + pairwise two-sided Mann-Whitney U computed on the data,
gene names in italic (fontstyle, so digits are italic too).

Gene list follows scanpy_CM_Fhl2.ipynb cells 43-48.

  ~/miniforge3/envs/scanpy/bin/python make_gene_violins_fixed.py
"""
import warnings, numpy as np, pandas as pd, scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from make_fig4_fixed import violin_panel, BASE, FIG

warnings.filterwarnings("ignore")
AD = f"{BASE}/adata_fixed/adata_Fhl2_analysed_fixed.h5ad"
OUT = f"{BASE}/adata_fixed"

GROUPS_OF_GENES = {
    "stress":       ["Nppa", "Nppb", "Xirp2", "Tlr4"],
    "myosin":       ["Mhrt", "Myh6", "Myh7"],
    "tf":           ["Mitf", "Nfia", "Nfe2l2", "Ppargc1a"],
    "mito_biogen":  ["Tfam", "Tfb1m", "Tfb2m", "Twnk", "Polrmt"],
    "mito_import":  ["Timm23", "Polg", "Tomm70a"],
    "oxphos":       ["Cox4i1", "Cox5b", "Cox6a1"],
}


def main():
    ad = sc.read_h5ad(AD)
    present = set(ad.raw.var_names)
    rows = []

    for tag, genes in GROUPS_OF_GENES.items():
        genes = [g for g in genes if g in present]
        missing = [g for g in GROUPS_OF_GENES[tag] if g not in present]
        if missing:
            print(f"  ({tag}: not in object -> {', '.join(missing)})")
        if not genes:
            continue
        # individual figures
        for g in genes:
            fig, ax = plt.subplots(figsize=(4, 4.6), dpi=200)
            rows.append(violin_panel(ax, ad, g, False, g, "expression", italic_title=True))
            fig.tight_layout()
            fig.savefig(f"{FIG}/violin_{g}_fixed.pdf", bbox_inches="tight")
            plt.close(fig)
        # one sheet per group
        fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4.6), dpi=200)
        axes = np.atleast_1d(axes)
        for ax, g in zip(axes, genes):
            violin_panel(ax, ad, g, False, g, "expression", italic_title=True)
        fig.tight_layout()
        fig.savefig(f"{FIG}/violins_{tag}_fixed.pdf", bbox_inches="tight")
        fig.savefig(f"{FIG}/violins_{tag}_fixed.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  {tag}: {', '.join(genes)}")

    stats = pd.DataFrame(rows).drop_duplicates(subset="feature")
    stats.to_csv(f"{OUT}/gene_violin_statistics_fixed.csv", index=False)
    pd.set_option("display.width", 160)
    print("\n--- computed statistics (per-cell KW + pairwise MWU) ---")
    print(stats.to_string(index=False))
    print("\nsaved to", FIG)


if __name__ == "__main__":
    main()
