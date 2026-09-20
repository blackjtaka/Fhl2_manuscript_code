#!/usr/bin/env python
"""Figure 6b — cardiomyocyte FHL2 by LVRR status, two-group model.

The legend states n = 4 per group and P = 0.035, so the panel is drawn from the
two-group model and nothing else. This reproduces branch [A] of
`recalc_fhl2_stats_260901.py`, which is where the reported P comes from:

  raw counts of CM_LVRR_DCM_120626.h5ad restricted to the barcodes retained in
  CM_analysed_130626.h5ad (8,447 cardiomyocyte nuclei from 8 donors)
  -> per-donor summed pseudobulk
  -> genes with >= 10 counts in >= 4 donors (15,820 genes)
  -> PyDESeq2, design ~LVRR with nonLVRR as reference
  -> FHL2 two-sided Wald test, nominal P reported (pre-specified single gene)

The earlier `Fig5_selfVAD_FHL2_recovery.py` drew a three-group model that also
included control donors C9/C10/C13/C15; there LVRR vs nonLVRR gives P = 0.017, which
is not the number in the legend. That script is kept in _excluded/superseded/.

Points are per-donor log1p(CP10K) of the donor-summed counts; the box shows the
median, the interquartile range, and whiskers to 1.5x the IQR, as the figure footer says.

  ~/miniforge3/envs/scanpy/bin/python make_fig6b_panel.py
"""
from pathlib import Path
import warnings, numpy as np, pandas as pd, anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy import sparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")

ADIR = Path("/Users/takahiro/Desktop/project/VAD_LVRRanalysis/scanpy/adata")
OUT = Path("/Users/takahiro/Desktop/reha_public_data/SelfVAD_LVRR/results")
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

LVRR = ["V33", "V64", "V70", "V84"]
NONLVRR = ["V30", "V51", "V56", "V109"]
C_LVRR, C_NON = "#2980b9", "#c0392b"

mpl.rcParams.update({
    "font.size": 7, "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def donor_sum(X, samples):
    donors, inv = np.unique(samples, return_inverse=True)
    ind = sparse.csr_matrix((np.ones(inv.shape[0]), (inv, np.arange(inv.shape[0]))),
                            shape=(len(donors), inv.shape[0]))
    agg = ind @ X
    return donors, (agg.toarray() if sparse.issparse(agg) else np.asarray(agg))


def versions():
    """DESeq2's dispersion fit depends on the optimiser, so the exact P moves between
    environments even though the effect estimate does not. Record what produced it."""
    import importlib.metadata as m, sys
    v = {"python": sys.version.split()[0]}
    for p in ("pydeseq2", "numpy", "scipy", "pandas", "anndata"):
        try:
            v[p] = m.version(p)
        except Exception:
            v[p] = "?"
    return v


def main():
    ver = versions()
    print("environment:", " | ".join(f"{k} {v}" for k, v in ver.items()))
    raw = ad.read_h5ad(ADIR / "CM_LVRR_DCM_120626.h5ad")
    ana = ad.read_h5ad(ADIR / "CM_analysed_130626.h5ad", backed="r")
    common = ana.obs_names.intersection(raw.obs_names)
    ana.file.close()
    raw = raw[common].copy()
    raw.obs["sample"] = raw.obs["sample"].astype(str)
    assert set(raw.obs["sample"]) == set(LVRR + NONLVRR), set(raw.obs["sample"])
    print(f"cardiomyocyte nuclei: {raw.n_obs:,} from {raw.obs['sample'].nunique()} donors")
    print("per-donor nuclei:", dict(raw.obs["sample"].value_counts()))

    genes = np.asarray(raw.var_names)
    donors, agg = donor_sum(sparse.csr_matrix(raw.X), raw.obs["sample"].values)
    cmat = pd.DataFrame(np.rint(agg).astype(int), index=donors, columns=genes)

    # ---- two-group DESeq2, nonLVRR as reference ---------------------------
    keep = (cmat >= 10).sum(0) >= 4
    c2 = cmat.loc[:, keep]
    meta = pd.DataFrame(
        {"LVRR": pd.Categorical(["LVRR" if s in LVRR else "nonLVRR" for s in c2.index],
                                categories=["nonLVRR", "LVRR"])}, index=c2.index)
    dds = DeseqDataSet(counts=c2, metadata=meta, design_factors="LVRR", n_cpus=1, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["LVRR", "LVRR", "nonLVRR"], quiet=True)
    st.summary()
    r = st.results_df.loc["FHL2"]
    p, lfc = float(r["pvalue"]), float(r["log2FoldChange"])
    print(f"\ngenes tested: {c2.shape[1]:,}")
    print(f"FHL2  log2FC = {lfc:+.4f}  lfcSE = {r['lfcSE']:.4f}  "
          f"two-sided Wald P = {p:.6g}  (BH across all genes: {r['padj']:.4g})")

    # ---- per-donor values for plotting ------------------------------------
    lib = agg.sum(1, keepdims=True); lib[lib == 0] = 1
    fi = int(np.where(genes == "FHL2")[0][0])
    pdm = pd.DataFrame({"donor": donors,
                        "FHL2": np.log1p(agg[:, fi] / lib.ravel() * 1e4)})
    pdm["group"] = np.where(pdm.donor.isin(LVRR), "LVRR", "non-LVRR")
    print("\nper-donor FHL2 (log1p CP10K):")
    print(pdm.sort_values(["group", "donor"]).to_string(index=False))

    # ---- panel -------------------------------------------------------------
    order, cols = ["LVRR", "non-LVRR"], [C_LVRR, C_NON]
    data = [pdm.loc[pdm.group == g, "FHL2"].values for g in order]

    fig, ax = plt.subplots(figsize=(2.1, 3.0), dpi=300)
    bp = ax.boxplot(data, positions=range(2), widths=0.55, patch_artist=True, whis=1.5,
                    medianprops={"color": "black", "linewidth": 0.9},
                    boxprops={"linewidth": 0.5}, whiskerprops={"linewidth": 0.5},
                    capprops={"linewidth": 0.5}, flierprops={"marker": ""})
    for box, c in zip(bp["boxes"], cols):
        box.set_facecolor(c); box.set_alpha(0.55)
    rng = np.random.default_rng(3)
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=16, color=cols[i],
                   alpha=0.95, edgecolor="black", lw=0.35, zorder=3)

    lo = min(v.min() for v in data); hi = max(v.max() for v in data)
    span = hi - lo
    y = hi + span * 0.12
    ax.plot([0, 0, 1, 1], [y, y + span * 0.05, y + span * 0.05, y], color="black", lw=0.6)
    ax.text(0.5, y + span * 0.07, f"{stars(p)}   $P$ = {p:.3f}",
            ha="center", va="bottom", fontsize=7)
    ax.set_ylim(lo - span * 0.18, y + span * 0.34)
    ax.set_xticks(range(2))
    ax.set_xticklabels([f"{g}\nn = {len(d)}" for g, d in zip(order, data)], fontsize=7)
    ax.set_ylabel("FHL2 (per-donor pseudobulk,\nlog1p CP10K)", fontsize=7)

    fig.tight_layout()
    fig.savefig(FIG / "Fig6b_FHL2_by_LVRR.pdf", bbox_inches="tight")
    fig.savefig(FIG / "Fig6b_FHL2_by_LVRR.png", dpi=300, bbox_inches="tight")
    pdm.to_csv(OUT / "Fig6b_per_donor_FHL2.csv", index=False)
    pd.DataFrame([{"contrast": "LVRR vs non-LVRR", "model": "~LVRR (two group)",
                   "n_per_group": 4, "n_nuclei": int(raw.n_obs),
                   "genes_tested": int(c2.shape[1]),
                   "log2FC": lfc, "lfcSE": float(r["lfcSE"]), "wald_p": p,
                   "padj_BH_all_genes": float(r["padj"]), **ver}]
                 ).to_csv(OUT / "Fig6b_statistics.csv", index=False)
    print("\nsaved:", FIG / "Fig6b_FHL2_by_LVRR.pdf")


if __name__ == "__main__":
    main()
