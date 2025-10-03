#!/usr/bin/env python3

import urllib.request
from datetime import datetime
from pathlib import Path
import csv
import statistics
import os

URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"

ts = datetime.now().strftime("%Y-%m-%d")

def download(url, dest):
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        if getattr(r, "status", 200) != 200:
            raise RuntimeError(f"HTTP error {r.status} for {url}")
        while True:
            chunk = r.read(8192)
            if not chunk:
                break
            f.write(chunk)

def parse_float(val):
    if not val or val.strip().lower() in {"na", "nan", "null", ""}:
        return None
    try:
        return float(val.strip())
    except ValueError:
        return None

def summarize_by_taxid(rows):
    """Only include records with both GC% and Size data."""
    by_taxid = {}

    for rec in rows:
        taxid = rec.get("TaxID")
        if not taxid:
            continue

        gc = parse_float(rec.get("GC%"))
        size_mb = parse_float(rec.get("Size (Mb)"))
        
        # Skip if either value is missing
        if gc is None or size_mb is None:
            continue
            
        size_bp = int(size_mb * 1_000_000)

        if taxid not in by_taxid:
            by_taxid[taxid] = {"gc": [], "length": []}
        
        by_taxid[taxid]["gc"].append(gc)
        by_taxid[taxid]["length"].append(size_bp)

    # Generate summary stats
    summary = []
    for taxid, data in by_taxid.items():
        n = len(data["gc"])  # Same as len(data["length"]) since we require both
        if n == 0:
            continue
            
        record = {
            "taxid": taxid,
            "n": n,
            "gc_mean": statistics.mean(data["gc"]),
            "gc_stdev": statistics.stdev(data["gc"]) if n > 1 else None,
            "length_mean": statistics.mean(data["length"]),
            "length_stdev": statistics.stdev(data["length"]) if n > 1 else None
        }
        summary.append(record)

    return summary

def write_csv(rows, path):
    if not rows:
        with open(path, "w", newline="") as f:
            f.write("taxid,n,gc_mean,gc_stdev,length_mean,length_stdev\n")
        return

    fieldnames = ["taxid", "n", "gc_mean", "gc_stdev", "length_mean", "length_stdev"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def main():
    src = 'eukaryotes.txt'

    print(f"Downloading {URL} -> {src}")
    download(URL, src)

    print("Parsing...")
    with open(src, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Parsed {len(rows)} records")

    print("Summarizing by TaxID...")
    summary = summarize_by_taxid(rows)
    summary_path = f"{ts}_NCBI_eukaryotes_stats.txt"
    write_csv(summary, summary_path)
    print(f"Wrote {len(summary)} taxa with complete data: {summary_path}")

    os.remove(src)

if __name__ == "__main__":
    main()