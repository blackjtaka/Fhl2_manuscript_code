#!/usr/bin/env python
"""Figure 6c — FHL2 in three public human cardiac cohorts (MAGNet, Kuppe 2022, Reichart 2022).

Self-contained: public inputs -> the six reported FHL2 contrasts -> the three box plots.

Inputs (place under ROOT; set REHA_PUBLIC_ROOT to override the default path)
  MAGNet/Counts.csv, MAGNet/phenoData.csv   https://github.com/mpmorley/MAGNet  (GSE141910)
  MAGNet/counts/GSM*_<sample>.csv.gz        GSE141910_RAW.tar, per-sample voom-normalised values (plotted)
  Kuppe2022/Kuppe2022_AllsnRNA.h5ad         CellxGene collection 8191c283-0816-424b-9b61-c3e1d6258a77
  Reichart2022/Reichart2022_CM.h5ad         CellxGene collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb

Statistics (FHL2 = ENSG00000115641, pre-specified single gene; nominal two-sided Wald P is reported,
the transcriptome-wide BH-adjusted value is written alongside it)
  MAGNet    raw counts, NF/DCM/HCM (PPCM excluded), samples with complete covariates, genes with >= 10 counts
            in >= 15 samples; PyDESeq2 ~ Library_Pool + RIN + age + gender + race + etiology, reference NF.
  Kuppe     cardiomyocytes (cell_type_original == 'Cardiomyocyte'); raw counts summed per
            donor x zone (major_labl) x patient_group; units with < 20 nuclei dropped;
            PyDESeq2 ~ patient_group, reference myogenic.
  Reichart  cardiomyocyte object; raw counts summed per donor; donors with < 20 nuclei dropped;
            PyDESeq2 ~ disease, reference normal.
  Kuppe/Reichart gene filter: >= 10 counts in >= max(3, n_units // 4) units.

Plotted values
  MAGNet    voom-normalised log2 FHL2 from the GEO per-sample files (all NF/DCM/HCM samples).
  Kuppe     mean of the object's log1p-normalised FHL2 over cardiomyocytes, per sample (>= 20 nuclei).
  Reichart  log1p(CP10K) of FHL2 in the donor-summed raw counts.

Run:  ~/miniforge3/envs/scanpy/bin/python fig6c_human_cohorts_FHL2.py
"""
import gzip  # noqa: F401  (pandas reads the .csv.gz files)
import importlib.metadata as md
import os
import re
import sys
import warnings
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("REHA_PUBLIC_ROOT", "/Users/takahiro/Desktop/reha_public_data"))
OUT = ROOT / "results" / "Fig6c"
FHL2_ENS = "ENSG00000115641"

C_NF = "#7d7d7d"; C_DCM = "#c0392b"; C_HCM = "#8e44ad"
C_MYO = "#7d7d7d"; C_FIB = "#d35400"; C_ISC = "#922b21"
C_NORM = "#7d7d7d"; C_ARVC = "#9b59b6"


def assert_integer_counts(X, what):
    data = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    assert data.size and np.all(data == np.round(data)) and data.min() >= 0, f"{what}: not non-negative integer counts"


def raw_counts(adata, what):
    assert adata.raw is not None, f"{what}: .raw (counts) is missing"
    assert_integer_counts(adata.raw.X, f"{what} .raw.X")
    return adata.raw.X, adata.raw.var.copy()


def fhl2_index(var, what):
    if FHL2_ENS in var.index:
        return int(np.where(var.index == FHL2_ENS)[0][0])
    assert "feature_name" in var.columns, f"{what}: no Ensembl index and no feature_name column"
    pos = np.where(var["feature_name"].astype(str).values == "FHL2")[0]
    assert len(pos) == 1, f"{what}: FHL2 not resolved uniquely"
    return int(pos[0])


def aggregate(adata, group_cols, what):
    X, var = raw_counts(adata, what)
    obs = adata.obs[group_cols].astype(str).copy()
    obs["__key__"] = obs[group_cols].agg("|".join, axis=1)
    keys, inv = np.unique(obs["__key__"].values, return_inverse=True)
    ind = sparse.csr_matrix((np.ones(adata.n_obs), (inv, np.arange(adata.n_obs))), shape=(len(keys), adata.n_obs))
    agg = (ind @ sparse.csr_matrix(X)).toarray()
    obs_df = pd.DataFrame([dict(zip(group_cols, k.split("|"))) for k in keys], index=keys)
    obs_df["n_cells"] = np.bincount(inv)
    return ad.AnnData(X=np.rint(agg).astype(int), obs=obs_df, var=var)


def run_de_FHL2(pdata, design, ref_lvl, contrasts):
    n = pdata.n_obs
    keep = (pdata.X >= 10).sum(axis=0) >= max(3, n // 4)
    sub = pdata[:, keep].copy()
    assert FHL2_ENS in sub.var_names, "FHL2 did not pass the expression filter"
    dds = DeseqDataSet(adata=sub, design=design, ref_level=ref_lvl, n_cpus=2, refit_cooks=True, quiet=True)
    dds.deseq2()
    out = []
    for contrast, name in contrasts:
        ds = DeseqStats(dds, contrast=contrast, n_cpus=2, quiet=True)
        ds.summary()
        r = ds.results_df.loc[FHL2_ENS]
        out.append({"contrast": name, "n_units": int(n), "genes_tested": int(sub.n_vars),
                    "log2FC": float(r["log2FoldChange"]), "lfcSE": float(r["lfcSE"]),
                    "wald_p_nominal": float(r["pvalue"]), "padj_BH_all_genes": float(r["padj"])})
    return out


# ---------------------------------------------------------------- MAGNet
def magnet():
    pheno = pd.read_csv(ROOT / "MAGNet/phenoData.csv").set_index("sample_name")
    rows = []
    for fp in sorted((ROOT / "MAGNet/counts").glob("*.csv.gz")):
        m = re.search(r"_([A-Z]\d+)\.csv\.gz$", fp.name)
        assert m, f"unexpected file name {fp.name}"
        df = pd.read_csv(fp, names=["gene", "value"], skiprows=1, header=None)
        v = df.loc[df["gene"] == FHL2_ENS, "value"]
        assert len(v) == 1, f"FHL2 not found once in {fp.name}"
        rows.append({"sample": m.group(1), "voom_FHL2": float(v.iloc[0])})
    per = pd.DataFrame(rows).merge(pheno[["etiology"]], left_on="sample", right_index=True, how="left")
    per = per[per["etiology"].isin(["NF", "DCM", "HCM"])].copy()

    counts = pd.read_csv(ROOT / "MAGNet/Counts.csv", index_col=0)
    samples = sorted(set(counts.columns) & set(pheno.index))
    counts, ph = counts[samples], pheno.loc[samples]
    keep = ph["etiology"].isin(["NF", "DCM", "HCM"])
    counts, ph = counts.loc[:, keep], ph.loc[keep]
    covar = ["etiology", "Library.Pool", "RIN", "age", "gender", "race"]
    ok = ph[covar].notna().all(axis=1)
    counts, ph = counts.loc[:, ok], ph.loc[ok]
    assert_integer_counts(counts.values, "MAGNet Counts.csv")
    counts = counts.loc[(counts >= 10).sum(axis=1) >= 15]
    print(f"MAGNet: {counts.shape[1]} samples x {counts.shape[0]} genes; {ph['etiology'].value_counts().to_dict()}")

    obs = ph[covar].copy()
    obs.columns = ["etiology", "Library_Pool", "RIN", "age", "gender", "race"]
    for c in ["etiology", "Library_Pool", "gender", "race"]:
        obs[c] = obs[c].astype(str)
    obs.index.name = "sample"
    adata = ad.AnnData(X=counts.T.values.astype(int), obs=obs, var=pd.DataFrame(index=[g.split(".")[0] for g in counts.index]))
    assert FHL2_ENS in adata.var_names
    dds = DeseqDataSet(adata=adata, design="~Library_Pool + RIN + age + C(gender) + C(race) + etiology",
                       ref_level=("etiology", "NF"), n_cpus=2, refit_cooks=True, quiet=True)
    dds.deseq2()
    res = []
    for case in ["DCM", "HCM"]:
        ds = DeseqStats(dds, contrast=["etiology", case, "NF"], n_cpus=2, quiet=True)
        ds.summary()
        r = ds.results_df.loc[FHL2_ENS]
        res.append({"contrast": f"MAGNet: {case}_vs_NF", "n_units": int(adata.n_obs), "genes_tested": int(adata.n_vars),
                    "log2FC": float(r["log2FoldChange"]), "lfcSE": float(r["lfcSE"]),
                    "wald_p_nominal": float(r["pvalue"]), "padj_BH_all_genes": float(r["padj"])})
    return per, res


# ---------------------------------------------------------------- Kuppe
def kuppe():
    a = sc.read_h5ad(ROOT / "Kuppe2022/Kuppe2022_AllsnRNA.h5ad")
    cm = a[a.obs["cell_type_original"].astype(str) == "Cardiomyocyte"].copy()
    del a
    print(f"Kuppe: {cm.n_obs:,} cardiomyocyte nuclei, {cm.obs['donor_id'].nunique()} donors, {cm.obs['sample'].nunique()} samples")
    # plotted values: mean of the object's log1p-normalised X per sample
    j = fhl2_index(cm.var, "Kuppe .var")
    x = np.asarray(cm.X[:, j].todense()).ravel() if sparse.issparse(cm.X) else np.asarray(cm.X[:, j]).ravel()
    assert x.max() < 20 and not np.all(x == np.round(x)), "Kuppe .X is expected to be log1p-normalised, not counts"
    o = cm.obs[["donor_id", "sample", "major_labl", "patient_group"]].astype(str).copy()
    o["FHL2_log1p"] = x
    per = (o.groupby(["donor_id", "sample", "major_labl", "patient_group"], observed=True)
             .agg(n_cells=("FHL2_log1p", "size"), mean_FHL2_log1p=("FHL2_log1p", "mean")).reset_index())
    per = per[per["n_cells"] >= 20].copy()

    pdata = aggregate(cm, ["donor_id", "major_labl", "patient_group"], "Kuppe")
    pdata = pdata[pdata.obs["n_cells"].astype(int) >= 20].copy()
    res = run_de_FHL2(pdata, "~patient_group", ("patient_group", "myogenic"),
                      [(["patient_group", "fibrotic", "myogenic"], "Kuppe: fibrotic_vs_myogenic"),
                       (["patient_group", "ischemic", "myogenic"], "Kuppe: ischemic_vs_myogenic")])
    return per, res


# ---------------------------------------------------------------- Reichart
def reichart():
    a = sc.read_h5ad(ROOT / "Reichart2022/Reichart2022_CM.h5ad")
    short = {"dilated cardiomyopathy": "DCM", "arrhythmogenic right ventricular cardiomyopathy": "ARVC",
             "non-compaction cardiomyopathy": "NCM", "normal": "Normal"}
    unknown = set(a.obs["disease"].astype(str)) - set(short)
    assert not unknown, f"Reichart: unmapped disease labels {unknown}"
    print(f"Reichart: {a.n_obs:,} cardiomyocyte nuclei, {a.obs['donor_id'].nunique()} donors")
    pdata = aggregate(a, ["donor_id", "disease"], "Reichart")
    del a
    pdata.obs["disease_short"] = pdata.obs["disease"].map(short)
    j = list(pdata.var_names).index(FHL2_ENS) if FHL2_ENS in pdata.var_names else fhl2_index(pdata.var, "Reichart .raw.var")
    lib = pdata.X.sum(axis=1)
    per = pd.DataFrame({"donor_id": pdata.obs["donor_id"].values, "disease": pdata.obs["disease_short"].values,
                        "n_cells": pdata.obs["n_cells"].values, "FHL2_log1p": np.log1p(pdata.X[:, j] / lib * 1e4)})
    assert per["donor_id"].is_unique, "Reichart: a donor carries more than one disease label"

    pdata = pdata[pdata.obs["n_cells"].astype(int) >= 20].copy()
    res = run_de_FHL2(pdata, "~disease_short", ("disease_short", "Normal"),
                      [(["disease_short", "DCM", "Normal"], "Reichart: DCM_vs_Normal"),
                       (["disease_short", "ARVC", "Normal"], "Reichart: ARVC_vs_Normal")])
    return per, res


# ---------------------------------------------------------------- figure
def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def box(ax, data, cols, point_size, seed):
    bp = ax.boxplot(data, positions=range(len(data)), widths=0.55, patch_artist=True,
                    medianprops={"color": "black", "linewidth": 0.8}, boxprops={"linewidth": 0.5},
                    whiskerprops={"linewidth": 0.5}, capprops={"linewidth": 0.5}, flierprops={"marker": "", "markersize": 0})
    for b, c in zip(bp["boxes"], cols):
        b.set_facecolor(c); b.set_alpha(0.55)
    rng = np.random.default_rng(seed)
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.13, 0.13, size=len(v)), v, s=point_size, color=cols[i],
                   alpha=0.85 if point_size > 5 else 0.45, edgecolor="black" if point_size > 5 else "none",
                   lw=0.3 if point_size > 5 else 0, zorder=3)


def panel(ax, data, labels, cols, p_list, ylabel, title, point_size, seed, pad):
    box(ax, data, cols, point_size, seed)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=6.3)
    ymax = max(v.max() for v in data); ymin = min(v.min() for v in data)
    for i, p in enumerate(p_list, start=1):
        ax.text(i, ymax + pad, stars(p), ha="center", fontsize=8)
    ax.set_ylim(ymin - 0.3, ymax + pad + 0.5)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold")


def main():
    for p in ["MAGNet/Counts.csv", "MAGNet/phenoData.csv", "MAGNet/counts", "Kuppe2022/Kuppe2022_AllsnRNA.h5ad", "Reichart2022/Reichart2022_CM.h5ad"]:
        assert (ROOT / p).exists(), f"missing input: {ROOT / p}"
    OUT.mkdir(parents=True, exist_ok=True)
    ver = {p: md.version(p) for p in ("pydeseq2", "numpy", "scipy", "pandas", "anndata", "scanpy")}
    ver["python"] = sys.version.split()[0]
    print("environment:", ver)

    mag_per, mag_res = magnet()
    kup_per, kup_res = kuppe()
    rei_per, rei_res = reichart()

    stats = pd.DataFrame(mag_res + kup_res + rei_res)
    stats = stats.assign(**ver)
    stats.to_csv(OUT / "Fig6c_FHL2_statistics.csv", index=False)
    mag_per.to_csv(OUT / "Fig6c_MAGNet_per_sample.csv", index=False)
    kup_per.to_csv(OUT / "Fig6c_Kuppe_per_sample.csv", index=False)
    rei_per.to_csv(OUT / "Fig6c_Reichart_per_donor.csv", index=False)
    print("\n" + stats[["contrast", "n_units", "genes_tested", "log2FC", "lfcSE", "wald_p_nominal",
                        "padj_BH_all_genes"]].to_string(index=False))

    P = stats.set_index("contrast")["wald_p_nominal"].to_dict()
    mpl.rcParams.update({"font.size": 7, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.0), dpi=300, gridspec_kw={"wspace": 0.45})
    o1 = ["NF", "DCM", "HCM"]
    d1 = [mag_per.loc[mag_per["etiology"] == e, "voom_FHL2"].values for e in o1]
    panel(axes[0], d1, [f"{e}\nn={len(v)}" for e, v in zip(o1, d1)], [C_NF, C_DCM, C_HCM],
          [P["MAGNet: DCM_vs_NF"], P["MAGNet: HCM_vs_NF"]], "FHL2 (voom log2)", "MAGNet (bulk RNA-seq)", 2, 0, 0.1)
    o2 = ["myogenic", "fibrotic", "ischemic"]
    d2 = [kup_per.loc[kup_per["patient_group"] == g, "mean_FHL2_log1p"].values for g in o2]
    n2 = [kup_per.loc[kup_per["patient_group"] == g, "donor_id"].nunique() for g in o2]
    panel(axes[1], d2, [f"{g}\nn={n} donors" for g, n in zip(o2, n2)], [C_MYO, C_FIB, C_ISC],
          [P["Kuppe: fibrotic_vs_myogenic"], P["Kuppe: ischemic_vs_myogenic"]], "FHL2 (log1p, per-sample mean CM)", "Kuppe 2022 (snRNA-seq)", 10, 1, 0.2)
    o3 = ["Normal", "DCM", "ARVC"]
    d3 = [rei_per.loc[rei_per["disease"] == d, "FHL2_log1p"].values for d in o3]
    panel(axes[2], d3, [f"{d}\nn={len(v)} donors" for d, v in zip(o3, d3)], [C_NORM, C_DCM, C_ARVC],
          [P["Reichart: DCM_vs_Normal"], P["Reichart: ARVC_vs_Normal"]], "FHL2 (log1p CP10K, per-donor CM)", "Reichart 2022 (snRNA-seq)", 10, 2, 0.2)
    fig.savefig(OUT / "Fig6c_FHL2_human_cohorts.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Fig6c_FHL2_human_cohorts.png", dpi=300, bbox_inches="tight")
    print("\nsaved:", OUT)


if __name__ == "__main__":
    main()
