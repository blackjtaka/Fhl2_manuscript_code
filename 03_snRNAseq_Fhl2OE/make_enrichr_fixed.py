#!/usr/bin/env python
"""Enrichment plots for the corrected Fhl2-OE run, in the original notebook's format.

Reproduces the layout used in scanpy_CM_Fhl2.ipynb (cells 102-108):
  one horizontal bar chart per gene-set library, adjusted P < 0.05, top 10 terms,
  sorted by -log10(adjusted P-value), bars in darkred, y axis inverted, 4 x 3 inches,
  term names with the trailing "(GO:…)" identifier stripped.

  ~/miniforge3/envs/scanpy/bin/python make_enrichr_fixed.py
"""
import warnings, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp

warnings.filterwarnings("ignore")
BASE = "/Users/takahiro/Desktop/project/Reha/snRNAseq_scanpy/Fhl2OE"
OUT, FIG = f"{BASE}/adata_fixed", f"{BASE}/figures_fixed"
GENESETS = ['MSigDB_Hallmark_2020', 'GO_Biological_Process_2023',
            'KEGG_2019_Mouse', 'Reactome_2022']
METRIC = '-log10(adjusted P-value)'


def run(genes, tag):
    genes = [g.strip() for g in genes]
    enr = gp.enrichr(gene_list=genes, gene_sets=GENESETS, organism='mouse', outdir=None)
    enr.results['n_genes'] = [int(x.split('/')[0]) for x in enr.results['Overlap']]
    enr.results.to_csv(f"{OUT}/enrichr_full_{tag}.csv", index=False)
    for geneset in GENESETS:
        df = enr.results[enr.results['Gene_set'] == geneset]
        df = df[df['Adjusted P-value'] < 0.05].copy()
        if df.empty:
            print(f"  {tag} / {geneset}: no term below 0.05")
            continue
        df[METRIC] = -np.log10(df['Adjusted P-value'])
        df = df.sort_values(METRIC, ascending=False)
        n_rank = 10
        plt.rcParams['axes.grid'] = False
        fig = plt.figure(figsize=(4, 3), dpi=200)
        plt.barh(width=df[:n_rank][METRIC],
                 y=[x.split(' (')[0] for x in df[:n_rank]['Term']],
                 color='darkred')
        plt.gca().invert_yaxis()
        plt.xlabel(METRIC)
        plt.title(f'{geneset}')
        plt.margins(y=0.02)
        plt.tight_layout()
        fname = f"enrichr_{tag}_{geneset}"
        fig.savefig(f"{FIG}/{fname}.pdf", bbox_inches="tight")
        fig.savefig(f"{FIG}/{fname}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  {tag} / {geneset}: {len(df)} terms, top = {df.iloc[0]['Term'][:55]}")
    return enr.results


def main():
    deg = pd.read_csv(f"{OUT}/DEG_Fhl2OE_vs_Ex_fixed.csv")
    up_fhl2 = deg.loc[deg.logfc > 0.5, "gene"]
    up_ex = deg.loc[deg.logfc < -0.5, "gene"]
    print(f"up in Fhl2 OE: {len(up_fhl2)} | up in Ex: {len(up_ex)}")
    run(up_fhl2.tolist(), "up_in_Fhl2OE")
    run(up_ex.tolist(), "up_in_Ex")
    print("saved to", FIG)


if __name__ == "__main__":
    main()
