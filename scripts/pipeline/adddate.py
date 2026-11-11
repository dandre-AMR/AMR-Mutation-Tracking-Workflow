#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path

BASE = Path("/mnt/c/AMR_Project")

INFILE  = BASE / "KP_Compiled_With_Sequence.tsv"
OUTFILE = BASE / "KP_Compiled_With_Sequence_And_Date.tsv"

ASSEMBLY_FILES = [
    BASE / "assembly_summary_refseq.txt",
    BASE / "assembly_summary_genbank.txt",
]

def load_date_map():
    """
    Build {GCA_accession -> seq_rel_date}
    """
    date_map = {}

    for f in ASSEMBLY_FILES:
        if not f.exists():
            continue

        with open(f, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 15:
                    continue

                acc = parts[0].strip()       # GCA_XXXXXXXXX.Y
                date = parts[14].strip()    # seq_rel_date

                if acc and date:
                    date_map[acc] = date

    return date_map

def extract_assembly(row, header):
    """
    Pull assembly accession from row.
    Prefers 'Assembly' column if present, otherwise infers from filename-style fields.
    """
    if "Assembly" in header:
        asm = row[header.index("Assembly")]
        return asm.split("_ASM")[0]

    raise RuntimeError("No Assembly column found in input file")

def main():
    date_map = load_date_map()
    print(f"[INFO] Loaded {len(date_map)} assembly dates")

    with open(INFILE, "r", encoding="utf-8") as ih, \
         open(OUTFILE, "w", encoding="utf-8", newline="") as oh:

        rdr = csv.reader(ih, delimiter="\t")
        wtr = csv.writer(oh, delimiter="\t")

        header = next(rdr)
        wtr.writerow(header + ["Assembly_Date"])

        added = 0

        for row in rdr:
            asm = extract_assembly(row, header)
            date = date_map.get(asm, "")
            if date:
                added += 1
            wtr.writerow(row + [date])

    print(f"[DONE] Dates added for {added} rows")

if __name__ == "__main__":
    main()
