#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import math

def log(m): print(f"[INFO] {m}", flush=True)

def clean_aa(s: pd.Series) -> pd.Series:
    s = s.fillna("").astype(str).str.upper()
    s = s.str.replace(r"\s+", "", regex=True)
    return s

def shannon_H(counts: pd.Series) -> float:
    N = counts.sum()
    if N <= 0: return 0.0
    p = counts / N
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

def evenness_E(counts: pd.Series) -> float:
    H = shannon_H(counts)
    S = int((counts > 0).sum())
    if S <= 1: return 0.0
    return float(H / math.log(S))

def simpson_D(counts: pd.Series) -> float:
    N = counts.sum()
    if N <= 0: return 0.0
    p = counts / N
    return float((p * p).sum())  # dominance

def p_tilde_max_jeffreys(xmax: int, N: int) -> float:
    # Jeffreys prior Beta(0.5, 0.5)
    if N < 0: return 0.0
    return float((xmax + 0.5) / (N + 1.0))

def compute_row(gene: str, variants: pd.Series) -> dict:
    # variants = AA_Sequence strings for this gene
    N_total = int(variants.shape[0])
    # keep non-empty sequences
    v = variants[variants != ""]
    S_variants = int(v.nunique()) if N_total else 0

    counts = v.value_counts()
    xmax = int(counts.max()) if not counts.empty else 0
    p_max = float(xmax / N_total) if N_total > 0 else 0.0
    p_tilde_max = p_tilde_max_jeffreys(xmax, N_total) if N_total > 0 else 0.0

    D = simpson_D(counts) if not counts.empty else 0.0
    SI = float(1.0 - D)                       # Gini–Simpson
    E = evenness_E(counts) if not counts.empty else 0.0
    OneMinusE = float(1.0 - E)

    # Redundancy: R = (N - S)/N
    Redundancy_R = float((N_total - S_variants)/N_total) if N_total > 0 else 0.0

    # Change_count: number of additional variants beyond the dominant one
    Change_count = int(max(S_variants - 1, 0))

    # Paper’s refined SI* — three-component product:
    SI_star = float(p_tilde_max * OneMinusE * Redundancy_R)

    # Also output the 2-component variant in case you want exactly that line:
    SI_star_2comp = float(p_tilde_max * OneMinusE)

    # Legacy % unique
    SI_old = float(100.0 * S_variants / N_total) if N_total > 0 else 0.0

    # Convenience combo you used earlier:
    SI_SimpsonR = float(SI * Redundancy_R)

    return {
        "Gene": gene,
        "N_total": N_total,
        "S_variants": S_variants,
        "p_max": p_max,
        "Change_count": Change_count,
        "p_tilde_max": p_tilde_max,
        "Evenness_E": E,
        "OneMinusE": OneMinusE,
        "Redundancy_R": Redundancy_R,
        "SI": SI,
        "Simpson_D": D,
        "SI_SimpsonR": SI_SimpsonR,
        "SI_old": SI_old,
        "SI_star": SI_star,
        "SI_star_2comp": SI_star_2comp
    }

def main():
    ap = argparse.ArgumentParser(description="Build KP-style stability table with SI* (dominance × (1-E) × redundancy).")
    ap.add_argument("--infile", required=True, help="Input CodonsAA TSV (e.g., AB_Compiled_With_Sequence_And_Date_CodonsAA.tsv)")
    ap.add_argument("--outfile", help="Output Excel (default: <infile>_Stability_KPHeaders.xlsx)")
    ap.add_argument("--gene-col", default="Element_symbol", help="Gene column (default: Element_symbol)")
    ap.add_argument("--aa-col", default="AA_Sequence", help="AA sequence column (default: AA_Sequence)")
    args = ap.parse_args()

    infile = Path(args.infile)
    outfile = Path(args.outfile) if args.outfile else infile.with_name(infile.stem + "_Stability_KPHeaders.xlsx")

    log(f"Reading {infile}")
    df = pd.read_csv(infile, sep="\t", low_memory=False)
    if args.gene_col not in df.columns or args.aa_col not in df.columns:
        raise SystemExit(f"Required columns not found: {args.gene_col}, {args.aa_col}")

    df = df.copy()
    df[args.aa_col] = clean_aa(df[args.aa_col])
    df = df[df[args.gene_col].notna() & (df[args.gene_col].astype(str).str.strip() != "")]

    rows = []
    for gene, g in df.groupby(args.gene_col, dropna=False):
        rows.append(compute_row(str(gene), g[args.aa_col]))

    out = pd.DataFrame(rows, columns=[
        "Gene","N_total","S_variants","p_max","Change_count","p_tilde_max",
        "Evenness_E","OneMinusE","Redundancy_R","SI","Simpson_D",
        "SI_SimpsonR","SI_old","SI_star","SI_star_2comp"
    ])

    # sort by SI* (three-component) desc, then N_total
    out = out.sort_values(["SI_star","N_total"], ascending=[False, False], ignore_index=True)

    log(f"Writing {outfile}")
    with pd.ExcelWriter(outfile, engine="openpyxl") as xw:
        out.to_excel(xw, index=False, sheet_name="Stability")
    log("Done.")

if __name__ == "__main__":
    main()
