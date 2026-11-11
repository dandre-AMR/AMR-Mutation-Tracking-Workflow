#!/usr/bin/env python3
import argparse, os, re
from pathlib import Path
import pandas as pd

AA20 = list("ARNDCEQGHILKMFPSTWYV")
AASET = set(AA20)

def log(m): print(f"[INFO] {m}", flush=True)

def infer_family(symbol: str) -> str:
    """Examples:
       'blaOXA-66'   -> 'blaOXA'
       'parC_S84L'   -> 'parC'
       'msr(E)'      -> 'msr(E)' (kept)
    """
    s = (symbol or "").strip()
    if not s: return s
    if "_" in s:
        return s.split("_", 1)[0]
    m = re.match(r"^([A-Za-z\(\)'\+\-]+?)(?:-?\d.*)?$", s)
    return m.group(1) if m and m.group(1) else s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True,
                    help="Master TSV with at least: Element symbol, AA_Sequence, Assembly, [Date], [Year]")
    ap.add_argument("--outfile-xlsx", default=None,
                    help="Output Excel (with sheets Wide, Long, Variant_Counts). Default: results/*_Raw_PerAssembly.xlsx")
    ap.add_argument("--outfile-wide-tsv", default=None,
                    help="Output Wide TSV. Default: results/*_AA_PerAssembly_Wide.tsv")
    ap.add_argument("--outfile-long-tsv", default=None,
                    help="Output Long TSV. Default: results/*_AA_PerAssembly_Long.tsv")
    ap.add_argument("--chunksize", type=int, default=200_000,
                    help="Row chunksize for streaming")
    args = ap.parse_args()

    in_path = Path(args.infile)
    if not in_path.exists():
        raise SystemExit(f"[ERROR] Missing infile: {in_path}")

    base_dir = in_path.parent.parent / "AA_Composition_Data" / "results"
    base_dir.mkdir(parents=True, exist_ok=True)

    # defaults
    stem = in_path.stem.replace(".tsv", "")
    out_xlsx = Path(args.outfile_xlsx) if args.outfile_xlsx else base_dir / f"AB_AB_AA_Raw_PerAssembly.xlsx"
    out_wide_tsv = Path(args.outfile_wide_tsv) if args.outfile_wide_tsv else base_dir / f"AB_AB_AA_PerAssembly_Wide.tsv"
    out_long_tsv = Path(args.outfile_long_tsv) if args.outfile_long_tsv else base_dir / f"AB_AB_AA_PerAssembly_Long.tsv"

    # Required columns (case sensitive as in your file)
    need_cols = ["Element symbol", "AA_Sequence", "Assembly"]
    opt_cols = ["Date", "Year"]
    usecols = need_cols + opt_cols

    log(f"Reading {in_path} in chunks…")
    agg = None
    n_rows = 0
    for chunk in pd.read_csv(in_path, sep="\t", usecols=lambda c: c in usecols, chunksize=args.chunksize):
        n_rows += len(chunk)

        # clean up
        chunk = chunk.dropna(subset=["AA_Sequence", "Assembly"])
        # family & variant
        chunk["Variant"] = chunk["Element symbol"].astype(str)
        chunk["Family"] = chunk["Element symbol"].astype(str).map(infer_family)
        # year fallback from Date if Year missing
        if "Year" not in chunk or chunk["Year"].isna().all():
            chunk["Year"] = pd.NA
        if "Date" in chunk:
            # if Year is NA, try first 4 of Date
            need_year = chunk["Year"].isna()
            y = pd.to_datetime(chunk.loc[need_year, "Date"], errors="coerce")
            chunk.loc[need_year, "Year"] = y.dt.year.astype("Int64")

        # Count AA per row (vectorized with str.count)
        seq = chunk["AA_Sequence"].astype(str)
        chunk["__len__"] = seq.str.len()
        for aa in AA20:
            chunk[aa] = seq.str.count(aa)

        # aggregate per (Family, Variant, Assembly, Year)
        keys = ["Family", "Variant", "Element symbol", "Assembly", "Year"]
        g = chunk.groupby(keys, dropna=False)[AA20 + ["__len__"]].sum().reset_index()

        agg = g if agg is None else pd.concat([agg, g], ignore_index=True)

    if agg is None or agg.empty:
        raise SystemExit("[ERROR] No data after filtering.")
    log(f"Aggregated rows: {len(agg)}")

    # proportions
    for aa in AA20:
        agg[aa] = agg[aa] / agg["__len__"].replace(0, pd.NA)

    # WIDE: proportions by AA columns
    wide_cols = ["Family", "Variant", "Element symbol", "Assembly", "Year"] + AA20
    wide = agg[wide_cols].copy()

    # LONG: melt AA cols to (AA, Proportion)
    long = wide.melt(
        id_vars=["Family", "Variant", "Element symbol", "Assembly", "Year"],
        value_vars=AA20,
        var_name="AA", value_name="Proportion"
    )

    # Variant counts (replicates per variant = #assemblies)
    var_counts = wide.groupby(["Family", "Variant"])["Assembly"].nunique().reset_index(name="n_assemblies")

    # write
    long.to_csv(out_long_tsv, sep="\t", index=False)
    wide.to_csv(out_wide_tsv, sep="\t", index=False)
    with pd.ExcelWriter(out_xlsx) as xw:
        wide.to_excel(xw, index=False, sheet_name="Wide")
        long.to_excel(xw, index=False, sheet_name="Long")
        var_counts.to_excel(xw, index=False, sheet_name="Variant_Counts")

    log(f"Wrote Wide TSV: {out_wide_tsv}")
    log(f"Wrote Long TSV: {out_long_tsv}")
    log(f"Wrote Excel: {out_xlsx}")

if __name__ == "__main__":
    main()
