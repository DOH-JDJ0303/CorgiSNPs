#!/usr/bin/env python3

"""
subtyper.py
Authors:
    Jared Johnson, jared.johnson@doh.wa.gov
    Zack Mudge, ZMudge@cdc.gov
"""

import argparse
import csv
import os
import sys
import sourmash as sm
import screed
import re
import logging
from typing import List, Tuple, Any, Optional

# ----------------------------
# Logging
# ----------------------------
def setup_logging(log_level: str = 'INFO', log_file: str | None = None) -> logging.Logger:
    level = getattr(logging, log_level.upper())
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(level)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(level)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger

# ----------------------------
# Utils
# ----------------------------
def sanitize_str(s: str) -> str:
    """Keep letters, numbers, underscores, dashes, and dots; replace others with '_'."""
    s = s.strip().replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9._-]', "_", s)
    s = re.sub(r'__+', "_", s)
    return s

def derive_taxon_from_gambit(path: str) -> Optional[str]:
    """
    Return predicted.name if predicted.rank == 'species' (case-insensitive).
    Otherwise return None.
    Uses the first non-empty row that has the required fields.
    """
    required_cols = {'predicted.name', 'predicted.rank'}
    with open(path, newline='', encoding='utf-8') as f:
        rdr = csv.DictReader(f)
        missing = required_cols - set(rdr.fieldnames or [])
        if missing:
            logging.warning(f"GAMBIT CSV missing required columns: {missing}")
            return None
        for row in rdr:
            name = (row.get('predicted.name') or '').strip()
            rank = (row.get('predicted.rank') or '').strip().lower()
            if not name:
                continue
            if rank == 'species':
                logging.info(f"GAMBIT predicted species: {name}")
                return name
            else:
                logging.info(f"GAMBIT top rank is not species (rank={row.get('predicted.rank')}); "
                             f"will return undefined subtype.")
                return None  # only consider the top/first row
    logging.warning("GAMBIT CSV had no usable rows.")
    return None

# ----------------------------
# I/O
# ----------------------------
def load_db_sig(sig_path: str) -> Tuple[List[Any], int, int]:
    """Load all signatures; ensure single ksize and scaled across the file."""
    sigs = list(sm.load_file_as_signatures(sig_path))
    if not sigs:
        raise ValueError(f"No signatures found in {sig_path}")

    kset = {s.minhash.ksize for s in sigs}
    sset = {s.minhash.scaled for s in sigs}
    if len(kset) > 1:
        raise ValueError(f"Multiple ksize values in {sig_path}: {kset}")
    if len(sset) > 1:
        raise ValueError(f"Multiple scaled values in {sig_path}: {sset}")

    logging.info(f"DB signatures: {len(sigs)}; ksize={list(kset)[0]}; scaled={list(sset)[0]}")
    return sigs, next(iter(kset)), next(iter(sset))

def load_seq(filepath: str, ksize: int = 31, scaled: int = 1000) -> sm.MinHash:
    """Build a MinHash from a sequence file (FASTA/FASTQ) by concatenating sequences."""
    seq = ''
    contigs = 0
    for rec in screed.open(filepath):
        seq += rec.sequence
        contigs += 1
    logging.info(f'Loaded {os.path.basename(filepath)}: {contigs:,} contigs; {len(seq):,} bp')

    mh = sm.MinHash(n=0, ksize=ksize, scaled=scaled)
    mh.add_sequence(seq, force=True)
    logging.info('Created MinHash from input sequence')
    return mh

def write_csv(path: str, row: dict) -> None:
    """Write a single-row CSV (overwrite)."""
    fields = [
        'sample', 'taxon', 'db', 'signature_file', 'ksize', 'scaled', 'n_db_signatures',
        'subtype', 'closest_subtype', 'closest_ani', 'threshold', 'passed_threshold'
    ]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerow(row)

# ----------------------------
# Main
# ----------------------------
def main():
    version = "v2.1.0"
    
    p = argparse.ArgumentParser(description='Subtyper for the pre-MycoSNP workflow')
    p.add_argument('--sample', required=True, help='Sample name')

    # taxon may be provided directly or derived from GAMBIT
    p.add_argument('--taxon', help='Taxon name (overrides --gambit if both provided)')
    p.add_argument('--gambit', help='GAMBIT output CSV (use predicted.name if predicted.rank=="species")')

    p.add_argument('--db', required=True, help='Path to the subtype database')
    p.add_argument('--seq', required=True, help='Path to the assembly or reads (FASTA/FASTQ)')
    p.add_argument('--distance_threshold', type=float, default=99.7,
                   help='Percent average nucleotide identity (ANI) threshold for a call')
    p.add_argument('--out', help='Output CSV path (default: <sample>.subtype.csv)')
    p.add_argument("--version", action="version", version=version)
    p.add_argument("--log-level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                   default='INFO', help="Logging level")
    p.add_argument("--log-file", help="Log file path (optional)")
    args = p.parse_args()

    setup_logging(args.log_level, args.log_file)

    sample = args.sample
    out_csv = args.out or f"{sample}_subtype.csv"
    db = args.db.rstrip('/')
    seq = args.seq

    # Resolve taxon (CLI > GAMBIT > None)
    taxon_val: Optional[str] = args.taxon
    if not taxon_val and args.gambit:
        logging.info(f"Deriving taxon from GAMBIT CSV: {args.gambit}")
        taxon_val = derive_taxon_from_gambit(args.gambit)

    # If we still don't have a species-level taxon, emit undefined and exit
    if not taxon_val:
        logging.warning("No species-level taxon available; writing undefined subtype.")
        write_csv(out_csv, {
            'sample': sample, 'taxon': args.taxon or '',
            'db': db, 'signature_file': '',
            'ksize': '', 'scaled': '', 'n_db_signatures': 0,
            'subtype': 'undefined', 'closest_subtype': '', 'closest_ani': '',
            'threshold': args.distance_threshold, 'passed_threshold': False
        })
        logging.info(f"Wrote: {out_csv}")
        return

    taxon_clean = sanitize_str(taxon_val)

    # manifest -> signature filepath
    manifest_path = os.path.join(db, "manifest.csv")
    if not os.path.isfile(manifest_path):
        raise ValueError(f"{manifest_path} does not exist")

    signature_rel = None
    with open(manifest_path, 'r', newline='', encoding='utf-8') as mf:
        reader = csv.DictReader(mf)
        for row in reader:
            if sanitize_str(row.get('taxon', '')) == taxon_clean:
                signature_rel = row.get('signature_filepath')
                break

    if not signature_rel:
        logging.warning(f"Taxon '{taxon_val}' not found in database manifest: {manifest_path}")
        write_csv(out_csv, {
            'sample': sample, 'taxon': taxon_val, 'db': db, 'signature_file': '',
            'ksize': '', 'scaled': '', 'n_db_signatures': 0,
            'subtype': 'undefined', 'closest_subtype': '', 'closest_ani': '',
            'threshold': args.distance_threshold, 'passed_threshold': False
        })
        logging.info(f"Wrote: {out_csv}")
        return

    signature_path = os.path.join(db, signature_rel)
    if not os.path.isfile(signature_path):
        raise ValueError(f"Signature file does not exist: {signature_path}")

    logging.info(f"Using signature file: {signature_path}")
    db_sigs, db_ksize, db_scaled = load_db_sig(signature_path)
    sample_mh = load_seq(seq, db_ksize, db_scaled)

    # Compare each DB signature to the single sample sketch
    top_ani = 0.0
    top_matches: set[str] = set()

    for sig in db_sigs:
        res = sig.minhash.containment_ani(sample_mh)
        try:
            dist = res.dist
        except AttributeError:
            dist = float(res)
        ani = 100.0 * (1.0 - dist)

        if ani > top_ani:
            top_ani = ani
            top_matches = {sig.name or "Undefined"}
        elif ani == top_ani:
            top_matches.add(sig.name or "Undefined")

    closest_subtype = " / ".join(sorted(top_matches)) if top_matches else ""
    passed = (top_ani >= float(args.distance_threshold)) and bool(top_matches)

    if not top_matches or not passed:
        subtype = 'undefined'
    else:
        subtype = closest_subtype

    row = {
        'sample': sample,
        'taxon': taxon_val,
        'db': db,
        'signature_file': signature_path,
        'ksize': db_ksize,
        'scaled': db_scaled,
        'n_db_signatures': len(db_sigs),
        'subtype': subtype,
        'closest_subtype': closest_subtype,
        'closest_ani': round(top_ani, 2) if top_matches else '',
        'threshold': args.distance_threshold,
        'passed_threshold': passed
    }
    write_csv(out_csv, row)
    logging.info(f"Wrote: {out_csv}")

if __name__ == '__main__':
    main()
