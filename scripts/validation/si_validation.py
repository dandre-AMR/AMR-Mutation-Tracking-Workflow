#!/usr/bin/env python3
"""
Stability Index (SI) validation and comparative figures.

Reads the per-gene stability tables written by scripts/stability/kp_stability_by_total.py
(one for K. pneumoniae, one for A. baumannii) and produces:

  1. Spearman correlations between SI and established diversity metrics (Table S4)
  2. A subsampling analysis of SI sensitivity to sample size (Figure S15)
  3. SI behaviour for two-variant genes (Figure S16) and SI vs Simpson (Figure S17)
  4. Variant proportion and variant-count panels (Figures 2 and 3)
  5. Shared-gene and shared-allele comparisons (Figures 5 and 6, Table S5)

SI is the two-component index used in the manuscript: SI = p~max x (1 - E),
column SI_star_2comp in the stability tables.

Example:
  python si_validation.py \
      --kp data/stability_tables/KP_Stability.xlsx \
      --ab data/stability_tables/AB_Stability.xlsx \
      --outdir results/validation
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Metal, biocide and heat-resistance genes and virulence genes are excluded from
# the comparative analyses (they remain in the stability tables).
EXCLUDE = re.compile(
    r"^(sil|pco|ars|mer|ter|cus|chr|czc|ncr|fieF|qac|clpK|hsp20|trxLHR|kefB|yfdX|shsP|psi|"
    r"ybt|clb|rmp|iuc|iut|iro|irp|fyu|mrk|peg|emrD)", re.I)
MIN_N = 50
SEED = 20261111
KP_COL, AB_COL, INK = "#2f6f9f", "#c0393b", "#222222"


def load(path):
    d = pd.read_excel(path)
    d["Gene"] = d["Gene"].astype(str)
    d["VP"] = d["S_variants"] / d["N_total"]
    d["SI2"] = d["SI_star_2comp"].clip(lower=0)
    d["SI3"] = d["SI_star"]
    d["AMR"] = ~d["Gene"].str.match(EXCLUDE)
    S = d["S_variants"].where(d["S_variants"] > 1, 1)
    d["H"] = d["Evenness_E"] * np.log(S)
    d["Hill1"] = np.exp(d["H"])
    d["Hill2"] = 1 / d["Simpson_D"]
    d["GiniSimpson"] = 1 - d["Simpson_D"]
    return d


def top_by_S(d, k):
    x = d[d["AMR"] & (d["N_total"] >= MIN_N)]
    return x.sort_values("S_variants", ascending=False).head(k)


def si_from_counts(counts):
    counts = counts[counts > 0]
    N, S = counts.sum(), len(counts)
    p = counts / N
    pt = (counts.max() + 0.5) / (N + 1)
    E = 0.0 if S == 1 else float(-(p * np.log(p)).sum() / np.log(S))
    return dict(SI2=pt * (1 - E), SI3=pt * (1 - E) * (N - S) / N, D=float((p ** 2).sum()), VP=S / N)


def synthetic(pmax, S, E):
    """Dominant variant at pmax, remaining variants Zipf-distributed with exponent fitted to E."""
    if S == 1:
        return np.array([1.0])
    def dist(a):
        w = np.arange(1, S) ** (-a)
        return np.concatenate([[pmax], w / w.sum() * (1 - pmax)])
    def f(a):
        p = dist(a)
        return -(p * np.log(p)).sum() / np.log(S) - E
    try:
        a = brentq(f, 0.0, 12)
    except ValueError:
        a = 0.0 if f(0) < 0 else 12
    return dist(a)


def correlations(kp, ab):
    cols = ["SI3", "p_max", "Simpson_D", "GiniSimpson", "H", "Hill1", "Hill2", "Evenness_E", "VP", "S_variants", "N_total"]
    out = {}
    for name, d in [("KP", kp), ("AB", ab)]:
        x = d[(d["N_total"] >= MIN_N) & d["AMR"]]
        out[name] = {"n": int(len(x))}
        for c in cols:
            r, p = spearmanr(x["SI2"], x[c])
            out[name][c] = {"rho": round(float(r), 3), "p": float(p)}
        out[name]["median_abs_SI2_minus_SI3"] = round(float((x["SI2"] - x["SI3"]).abs().median()), 4)
    return out


def subsampling(kp, ab, outdir, reps=300):
    rng = np.random.default_rng(SEED)
    genes = [("KP", "sul2"), ("KP", "aph(6)-Id"), ("KP", "blaOXA"), ("KP", "oqxB"), ("KP", "blaTEM"),
             ("AB", "ant(3'')-IIa"), ("AB", "blaADC"), ("AB", "blaTEM-1")]
    Ns = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000]
    rows = []
    for sp, g in genes:
        d = kp if sp == "KP" else ab
        hit = d[d["Gene"] == g]
        if hit.empty:
            continue
        r = hit.iloc[0]
        p = synthetic(r["p_max"], int(r["S_variants"]), r["Evenness_E"])
        Ntot = int(r["N_total"])
        pop = np.floor(p * Ntot).astype(int)
        pop[0] += Ntot - pop.sum()
        full = si_from_counts(pop.astype(float))
        for N in Ns:
            if N > Ntot:
                continue
            sims = [si_from_counts(rng.multivariate_hypergeometric(pop, N).astype(float)) for _ in range(reps)]
            for m in ["SI2", "SI3", "D", "VP"]:
                v = np.array([s[m] for s in sims])
                rows.append(dict(sp=sp, gene=g, N=N, metric=m, full=full[m], mean=v.mean(),
                                 lo=np.percentile(v, 2.5), hi=np.percentile(v, 97.5)))
    sim = pd.DataFrame(rows)
    sim.to_csv(outdir / "subsampling.csv", index=False)
    s = sim[sim.metric == "SI2"]
    fig, axs = plt.subplots(2, 4, figsize=(12, 6), sharex=True, sharey=True)
    for ax, (key, q) in zip(axs.flat, s.groupby(["sp", "gene"], sort=False)):
        c = KP_COL if key[0] == "KP" else AB_COL
        ax.fill_between(q.N, q.lo, q.hi, color=c, alpha=0.18, lw=0)
        ax.plot(q.N, q["mean"], color=c, lw=2, marker="o", ms=4)
        ax.axhline(q.full.iloc[0], color="#444", ls="--", lw=1)
        ax.set_xscale("log"); ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{key[1]} ({'Kp' if key[0] == 'KP' else 'Ab'})", fontstyle="italic", loc="left", fontsize=10)
    for ax in axs[1]:
        ax.set_xlabel("Subsample size (N)")
    for ax in axs[:, 0]:
        ax.set_ylabel("SI")
    fig.tight_layout(); fig.savefig(outdir / "FigureS15_subsampling.png", dpi=300); plt.close(fig)
    summary = (s.assign(bias=s["mean"] - s["full"]).groupby("N")["bias"]
               .agg(mean_bias="mean", max_abs_bias=lambda b: b.abs().max()).round(3))
    return summary.reset_index().to_dict("records")


def two_variant_figure(ab, outdir):
    f = np.linspace(0.5, 0.999, 200)
    E = -(f * np.log(f) + (1 - f) * np.log(1 - f)) / np.log(2)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(f, f * (1 - E), color=INK, lw=2, label="SI")
    ax.plot(f, f ** 2 + (1 - f) ** 2, color=KP_COL, lw=2, ls="--", label="Simpson concentration")
    ax.plot(f, f, color="#777", lw=1.5, ls=":", label="Berger-Parker dominance")
    pts = ab[(ab["S_variants"] == 2) & (ab["N_total"] >= MIN_N)]
    ax.scatter(pts["p_max"], pts["SI2"], s=22, color=AB_COL, zorder=3, label=f"A. baumannii, S = 2, N >= {MIN_N}")
    ax.set_xlabel("Frequency of the dominant variant (S = 2)"); ax.set_ylabel("Index value")
    ax.set_xlim(0.5, 1); ax.set_ylim(0, 1.02); ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig(outdir / "FigureS16_two_variants.png", dpi=300); plt.close(fig)


def si_vs_simpson_figure(kp, ab, outdir):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for ax, d, nm, c in [(axs[0], kp, "K. pneumoniae", KP_COL), (axs[1], ab, "A. baumannii", AB_COL)]:
        x = d[(d["N_total"] >= MIN_N) & d["AMR"]]
        for cond, mk, lab, col in [(x.S_variants == 1, "s", "S = 1", "#999"),
                                   (x.S_variants.between(2, 3), "^", "S = 2-3", c),
                                   (x.S_variants > 3, "o", "S > 3", "#333")]:
            q = x[cond]
            ax.scatter(q["Simpson_D"], q["SI2"], s=14, marker=mk, color=col, alpha=0.75, label=f"{lab} (n = {len(q)})")
        ax.set_title(nm, fontstyle="italic", loc="left"); ax.set_xlabel("Simpson concentration")
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    axs[0].set_ylabel("Stability Index (SI)")
    fig.tight_layout(); fig.savefig(outdir / "FigureS17_SI_vs_Simpson.png", dpi=300); plt.close(fig)


def species_figure(d, color, path):
    t = top_by_S(d, 20)
    a = t.sort_values("VP", ascending=False)
    fig, axs = plt.subplots(1, 2, figsize=(11, 6.2))
    for ax, lab, vals, xl, fmt, panel in [
            (axs[0], a["Gene"], 100 * a["VP"], "Variant proportion, VP = S/N (%)", "{:.1f}", "A"),
            (axs[1], t["Gene"], t["S_variants"], "Distinct sequence variants (S)", "{:,.0f}", "B")]:
        y = np.arange(len(lab))[::-1]
        ax.barh(y, vals, color=color, height=0.72)
        ax.set_yticks(y); ax.set_yticklabels(lab, fontstyle="italic")
        for yi, v in zip(y, vals):
            ax.text(v + vals.max() * 0.01, yi, fmt.format(v), va="center", fontsize=8.5)
        ax.set_xlim(0, vals.max() * 1.15); ax.set_xlabel(xl)
        ax.set_title(panel, loc="left", fontweight="bold")
    fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig)


def shared_genes(kp, ab, outdir):
    a50 = set(top_by_S(ab, 50)["Gene"])
    rows = []
    for g in top_by_S(kp, 50)["Gene"]:
        if g in a50:
            k = kp[kp.Gene == g].iloc[0]; a = ab[ab.Gene == g].iloc[0]
            rows.append(dict(gene=g, SI_KP=k.SI2, SI_AB=a.SI2, N_KP=int(k.N_total), N_AB=int(a.N_total),
                             S_KP=int(k.S_variants), S_AB=int(a.S_variants)))
    sh = pd.DataFrame(rows)
    sh["dSI_KP_minus_AB"] = sh.SI_KP - sh.SI_AB
    sh = sh.sort_values("dSI_KP_minus_AB", ascending=False)
    sh.to_csv(outdir / "shared_genes.csv", index=False)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 8), sharey=True)
    y = np.arange(len(sh))[::-1]
    axL.barh(y, sh.SI_KP, color=KP_COL); axL.set_xlim(1.05, 0); axL.set_title("A   K. pneumoniae", loc="left", fontstyle="italic")
    axR.barh(y, sh.SI_AB, color=AB_COL); axR.set_xlim(0, 1.05); axR.set_title("B   A. baumannii", loc="left", fontstyle="italic")
    for ax in (axL, axR):
        for x in (0.1, 0.5):
            ax.axvline(x, color="#aaa", ls="--", lw=0.8)
        ax.set_xlabel("Stability Index (SI)")
    axL.set_yticks(y); axL.set_yticklabels(sh.gene, fontstyle="italic")
    fig.tight_layout(); fig.savefig(outdir / "Figure5_shared_genes.png", dpi=300); plt.close(fig)
    return sh


def shared_alleles(kp, ab, outdir):
    abg = ab.set_index("Gene")
    rows = []
    for _, k in kp.iterrows():
        g = k.Gene
        if g in abg.index and g.startswith("bla") and re.search(r"-\d", g):
            a = abg.loc[g]
            if k.N_total >= MIN_N and a.N_total >= MIN_N:
                rows.append(dict(allele=g, N_KP=int(k.N_total), S_KP=int(k.S_variants), pmax_KP=k.p_max,
                                 SI_KP=k.SI2, D_KP=k.Simpson_D, N_AB=int(a.N_total), S_AB=int(a.S_variants),
                                 pmax_AB=a.p_max, SI_AB=a.SI2, D_AB=a.Simpson_D))
    al = pd.DataFrame(rows).sort_values("N_KP", ascending=False)
    al["dSI_KP_minus_AB"] = al.SI_KP - al.SI_AB
    al.round(3).to_csv(outdir / "shared_alleles.csv", index=False)
    return al


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kp", required=True, help="K. pneumoniae stability table (.xlsx)")
    ap.add_argument("--ab", required=True, help="A. baumannii stability table (.xlsx)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--reps", type=int, default=300, help="Subsamples per sample size (default 300)")
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    kp, ab = load(args.kp), load(args.ab)

    res = {"correlations": correlations(kp, ab)}
    res["subsampling_SI_bias_by_N"] = subsampling(kp, ab, outdir, args.reps)
    two_variant_figure(ab, outdir)
    si_vs_simpson_figure(kp, ab, outdir)
    species_figure(kp, KP_COL, outdir / "Figure2_KP.png")
    species_figure(ab, AB_COL, outdir / "Figure3_AB.png")
    sh = shared_genes(kp, ab, outdir)
    al = shared_alleles(kp, ab, outdir)
    res["n_shared_genes"] = int(len(sh))
    res["n_shared_genes_higher_in_KP"] = int((sh.dSI_KP_minus_AB > 0).sum())
    res["n_shared_alleles"] = int(len(al))
    for name, d in [("KP", kp), ("AB", ab)]:
        x = d[(d.N_total >= MIN_N) & d.AMR]
        res[f"{name}_genes_N_ge_{MIN_N}"] = int(len(x))
        res[f"{name}_single_variant"] = int((x.S_variants == 1).sum())
        res[f"{name}_two_variants"] = int((x.S_variants == 2).sum())
        res[f"{name}_two_variants_pmax_lt_0.6"] = int(((x.S_variants == 2) & (x.p_max < 0.6)).sum())
    (outdir / "summary.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "correlations"}, indent=1))


if __name__ == "__main__":
    main()
