#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch KP or AB genomes (FASTA) from NCBI RefSeq/GenBank with parallel workers.

Examples:
  python fetch_all_kp_ab.py /mnt/c/AMR_Project --species KP --threads 8
  python fetch_all_kp_ab.py /mnt/c/AMR_Project --species AB --source refseq --threads 4

Outputs under basedir:
  KP_FASTA/ or AB_FASTA/
  kp_accessions.txt or ab_accessions.txt
"""

import argparse
import csv
import gzip
import io
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Tuple
from urllib.request import urlopen, Request

REFSEQ_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
GENBANK_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"

# organism filters (case-insensitive)
ORG_PATTERNS = {
    "KP": re.compile(r"\bKlebsiella\s+pneumoniae\b", re.IGNORECASE),
    "AB": re.compile(r"\bAcinetobacter\s+baumannii\b", re.IGNORECASE),
}

def http_get(url: str, timeout: int = 60) -> bytes:
    """Download a URL and return raw bytes with a sensible UA and timeout."""
    req = Request(url, headers={"User-Agent": "ncbi-fetch/1.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()

def read_assembly_summary(url: str) -> List[dict]:
    """Fetch and parse an NCBI assembly_summary_*.txt as list of dicts (TSV)."""
    raw = http_get(url)
    text = raw.decode("utf-8", errors="replace").splitlines()
    rows = []
    reader = csv.reader(text, delimiter="\t")
    headers = None
    for line in reader:
        if not line:
            continue
        if line[0].startswith("#"):
            # header line is the one that starts with '#'
            if "# assembly_accession" in line[0] or any("assembly_accession" in x for x in line):
                # Normalize headers: strip leading '# ' if present
                headers = [h.lstrip("# ").strip() for h in line]
            continue
        if headers is None:
            # Skip lines until we saw header
            continue
        if len(line) < len(headers):
            # pad (rare)
            line += [""] * (len(headers) - len(line))
        row = {h: v for h, v in zip(headers, line)}
        rows.append(row)
    return rows

def to_https(ftp_path: str) -> str:
    """Convert NCBI ftp path to https path."""
    return ftp_path.replace("ftp://", "https://")

def genomic_fna_url_from_ftp(ftp_path: str) -> str:
    """
    Build the *_genomic.fna.gz URL from assembly ftp path.
    Example ftp path:
      ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../GCF_XXXX.1_asmdir
    """
    https = to_https(ftp_path.rstrip("/"))
    base = https.split("/")[-1]  # e.g., GCF_XXXX.1_asmdir
    return f"{https}/{base}_genomic.fna.gz"

def filter_species(rows: List[dict], species_code: str) -> List[dict]:
    pat = ORG_PATTERNS[species_code]
    keep = []
    for r in rows:
        org = r.get("organism_name", "")
        if pat.search(org or ""):
            keep.append(r)
    return keep

def latest_per_accession(rows: List[dict]) -> List[dict]:
    """
    Keep only rows where 'version_status' == 'latest'.
    """
    out = [r for r in rows if (r.get("version_status") or "").lower() == "latest"]
    return out

def prepare_tasks(rows: List[dict]) -> List[Tuple[str, str]]:
    """
    From assembly_summary rows, construct list of (accession, fna_url).
    """
    tasks = []
    for r in rows:
        acc = r.get("assembly_accession", "").strip()
        ftp_path = r.get("ftp_path", "").strip()
        if not acc or not ftp_path or ftp_path == "na":
            continue
        url = genomic_fna_url_from_ftp(ftp_path)
        tasks.append((acc, url))
    return tasks

def write_accessions(path: Path, tasks: List[Tuple[str, str]]):
    with path.open("w", encoding="utf-8") as fh:
        for acc, url in tasks:
            fh.write(f"{acc}\t{url}\n")

def download_one(url: str, out_fna_gz: Path, retries: int = 3, backoff: float = 2.0) -> bool:
    """
    Download a single .fna.gz to out_fna_gz.
    Returns True if file exists at end (downloaded or already present).
    """
    if out_fna_gz.exists() and out_fna_gz.stat().st_size > 0:
        return True
    attempt = 0
    while attempt <= retries:
        try:
            data = http_get(url, timeout=120)
            out_fna_gz.parent.mkdir(parents=True, exist_ok=True)
            with out_fna_gz.open("wb") as fh:
                fh.write(data)
            return True
        except Exception:
            attempt += 1
            if attempt > retries:
                return False
            time.sleep(backoff * attempt)
    return False

def gunzip_to_fna(src_gz: Path, dst_fna: Path) -> bool:
    """
    Decompress .gz to .fna (skip if exists and non-empty).
    """
    if dst_fna.exists() and dst_fna.stat().st_size > 0:
        return True
    try:
        with gzip.open(src_gz, "rb") as gzf, dst_fna.open("wb") as out:
            out.write(gzf.read())
        return True
    except Exception:
        return False

def worker_task(acc: str, url: str, outdir: Path) -> Tuple[str, bool, str]:
    """
    Worker: download acc to outdir, return (acc, success, message).
    """
    try:
        base = url.split("/")[-1]  # *_genomic.fna.gz
        out_gz = outdir / base
        out_fna = outdir / base.replace("_genomic.fna.gz", "_genomic.fna")
        ok = download_one(url, out_gz)
        if not ok:
            return (acc, False, "download_failed")
        if not gunzip_to_fna(out_gz, out_fna):
            return (acc, False, "gunzip_failed")
        return (acc, True, "ok")
    except Exception as e:
        return (acc, False, f"exception: {e}")

def fetch_species(basedir: Path, species_code: str, source: str, threads: int):
    """
    Download genomes for species_code: 'KP' or 'AB'
    source: 'refseq', 'genbank', or 'both'
    """
    species_label = {"KP": "Klebsiella pneumoniae", "AB": "Acinetobacter baumannii"}[species_code]
    fasta_dir = basedir / ("KP_FASTA" if species_code == "KP" else "AB_FASTA")
    acc_path = basedir / ("kp_accessions.txt" if species_code == "KP" else "ab_accessions.txt")

    # get assembly summaries
    rows: List[dict] = []
    if source in ("refseq", "both"):
        rows += read_assembly_summary(REFSEQ_SUMMARY)
    if source in ("genbank", "both"):
        rows += read_assembly_summary(GENBANK_SUMMARY)

    # filter and prep
    rows = filter_species(rows, species_code)
    rows = latest_per_accession(rows)
    tasks = prepare_tasks(rows)
    tasks = sorted(tasks, key=lambda x: x[0])  # stable order

    # write accessions list
    write_accessions(acc_path, tasks)

    # download in parallel
    print(f"[INFO] {species_label}: {len(tasks)} assemblies to download from {source.upper()}.")
    fasta_dir.mkdir(parents=True, exist_ok=True)

    ok, fail = 0, 0
    with ThreadPoolExecutor(max_workers=threads) as exe:
        futs = [exe.submit(worker_task, acc, url, fasta_dir) for acc, url in tasks]
        for fut in as_completed(futs):
            acc, success, msg = fut.result()
            if success:
                ok += 1
                if ok % 50 == 0:
                    print(f"[OK] {ok} done...")
            else:
                fail += 1
                print(f"[WARN] {acc}: {msg}")

    print(f"[DONE] {species_label}: success={ok}, failed={fail}, outdir={fasta_dir}")

def main():
    ap = argparse.ArgumentParser(description="Download KP or AB genomes from NCBI with parallel workers.")
    ap.add_argument("basedir", help="Base directory (e.g., /mnt/c/AMR_Project)")
    ap.add_argument("--species", choices=["KP", "AB"], required=True,
                    help="Which species to download (KP=Klebsiella pneumoniae, AB=Acinetobacter baumannii)")
    ap.add_argument("--source", choices=["refseq", "genbank", "both"], default="both",
                    help="Assembly source(s) to use")
    ap.add_argument("--threads", type=int, default=os.cpu_count(),
                    help="Number of parallel workers (default: all cores)")
    args = ap.parse_args()

    basedir = Path(args.basedir)
    basedir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Base: {basedir}")
    print(f"[INFO] Species: {args.species}")
    print(f"[INFO] Source: {args.source}")
    print(f"[INFO] Threads: {args.threads}")

    fetch_species(basedir, args.species, args.source, args.threads)

if __name__ == "__main__":
    main()
