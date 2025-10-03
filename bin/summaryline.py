#!/usr/bin/env python3

import json
import csv
import screed
import argparse
import logging
from typing import List, Dict, Any, Optional


# ----------------------------
# File helper functions
# ----------------------------

def load_json(path: str, source: Optional[str] = None) -> Dict[str, Any]:
    """Load JSON file and return as dict."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    logging.info(f"Loaded JSON with {len(data)} top-level keys ({source})")
    return data


def load_csv(path: str, source: Optional[str] = None) -> List[Dict[str, str]]:
    """Load CSV into a list of dictionaries (rows)."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        data = [row for row in reader]
    logging.info(f"Loaded CSV with {len(data)} rows ({source})")
    return data


def load_fasta(path: str):
    """Return a context manager over screed records for streaming."""
    return screed.open(path)


def z_score(value: float, mean: float, sd: float) -> Optional[float]:
    # z-score is valid for any mean; only sd must be > 0
    if sd and sd > 0:
        return (value - mean) / sd
    return None


def genome_stats(records) -> Dict[str, Any]:
    """
    Stream through FASTA once; avoid building a giant string and avoid .get on screed records.
    """
    contigs = 0
    length = 0
    gc_count = 0
    for rec in records:
        seq = rec.sequence.upper()
        contigs += 1
        length += len(seq)
        gc_count += (seq.count('G') + seq.count('C'))
    gc = (100.0 * gc_count / length) if length > 0 else None
    logging.info(f"Loaded FASTA with {contigs} contigs")
    return {'denovo_contigs': contigs, 'denovo_length': length, 'denovo_gc': gc}


def taxid_stats(data: Dict[str, Any],
                ncbi_stats: List[Dict[str, str]],
                taxid: Any,
                min_n: int) -> Dict[str, Any]:
    """
    Attach z-scores for assembly length and GC if ncbi_stats has enough samples.
    Prefer exact taxid match; fall back to case-insensitive name match.
    """
    res: Dict[str, Any] = {}

    name = (data.get('species') or '').strip()
    length = data.get('denovo_length', 0)
    gc = data.get('denovo_gc', None)
    if not length or gc is None:
        return res

    wanted_taxid = (str(taxid).strip() if taxid is not None else None)

    for rec in ncbi_stats:
        rec_taxid = (str(rec.get('taxid') or '').strip() or None)
        rec_name  = (rec.get('name') or '').strip()

        taxid_match = bool(wanted_taxid and rec_taxid and rec_taxid == wanted_taxid)
        name_match  = bool(name and rec_name and rec_name.lower() == name.lower())
        if not (taxid_match or name_match):
            continue

        length_mean = float(rec.get("length_mean") or 0)
        length_sd   = float(rec.get("length_stdev") or 0)
        gc_mean     = float(rec.get("gc_mean") or 0)
        gc_sd       = float(rec.get("gc_stdev") or 0)
        n           = int(float(rec.get("n") or 0))

        if n >= min_n:
            res['denovo_length_z'] = z_score(float(length), length_mean, length_sd)
            res['denovo_gc_z']     = z_score(float(gc), gc_mean, gc_sd)
        break

    return res


def parse_read_stats(stats_dict: Dict[str, Any]) -> Dict[str, Any]:
    res: Dict[str, Any] = {}
    for k, v in stats_dict.get('summary', {}).items():
        if k not in ['before_filtering', 'after_filtering']:
            continue
        for k2, v2 in v.items():
            res[f"{k2}_{k}"] = v2
    return res


def parse_species(species_dict: Dict[str, str]):
    name = species_dict.get('predicted.name')
    rank = species_dict.get('predicted.rank')
    taxid = species_dict.get('predicted.ncbi_id')

    conf = 'high' if rank == 'species' else 'low'

    if not name:
        name = species_dict.get('next.name')

    res = {'species': name, 'species_confidence': conf}
    return res, taxid


def parse_subtype(subtype_dict: Dict[str, str]) -> Dict[str, Any]:
    name = subtype_dict.get('subtype')
    conf = subtype_dict.get('closest_ani')
    if name:
        name = name.split('-', 1)[0]
    return {'subtype': name, 'subtype_ani': conf}


def auto_qc(data: Dict[str, Any],
            min_depth: int,
            min_qual: float,
            max_z: float) -> Dict[str, Any]:

    qc_status = 'PASS'
    qc_reason: List[str] = []

    # Require only the essentials; z-scores are optional below
    req_cols = ['species', 'subtype', 'denovo_depth', 'q30_rate_after_filtering']
    for col in req_cols:
        if col not in data:
            return data | {'qc_status': 'UNKNOWN',
                           'qc_reason': [f'{col} not determined']}

    # Quality checks
    try:
        if float(data['q30_rate_after_filtering']) < float(min_qual):
            qc_status = 'FAIL'
            qc_reason.append(f'Q30 Rate < {min_qual}')
    except (TypeError, ValueError):
        qc_status = 'UNKNOWN'
        qc_reason.append('Invalid Q30 rate')

    try:
        if int(float(data['denovo_depth'])) < int(min_depth):
            qc_status = 'FAIL'
            qc_reason.append(f'Assembly Depth < {min_depth}')
    except (TypeError, ValueError):
        qc_status = 'UNKNOWN'
        qc_reason.append('Invalid assembly depth')

    # Optional z-score checks
    for key, label in [('denovo_length_z', 'Unusual assembly length'),
                       ('denovo_gc_z', 'Unusual GC content')]:
        val = data.get(key)
        try:
            if val is not None and abs(float(val)) >= float(max_z):
                qc_status = 'FAIL'
                qc_reason.append(label)
        except (TypeError, ValueError):
            qc_status = 'UNKNOWN'
            qc_reason.append(f'Invalid {key}')

    return data | {'qc_status': qc_status, 'qc_reason': qc_reason}


def parse_samplesheet(data: List[Dict[str, str]], sample: str) -> Dict[str, Any]:
    found = False
    out: Dict[str, Any] = {}
    for rec in data:
        rec_name = rec.get('sample')
        if not rec_name or rec_name != sample:
            continue
        found = True
        if rec.get('species'):
            out |= {'species': rec['species'], 'species_confidence': 'Manual'}
        if rec.get('subtype'):
            out |= {'subtype': rec['subtype']}
        if rec.get('reference'):
            out |= {'reference': rec['reference'], 'reference_source': 'Manual'}
    if not found:
        logging.warning(f"{sample} not found in samplesheet!")
    elif out:
        logging.info(f"Using values from samplesheet: {list(out.keys())}")
    else:
        logging.info(f"No relevant values found in samplesheet for {sample}")
    return out


# ----------------------------
# Main
# ----------------------------

def main():
    version = "1.0"

    parser = argparse.ArgumentParser(
        description="Summarize outputs from various workflows"
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--samplesheet", required=True, help="Samplesheet path")
    parser.add_argument("--ncbi_stats", default=None, help="NCBI stats file.")
    parser.add_argument("--read_stats", help="Fastp summary output (JSON)")
    parser.add_argument("--species", help="GAMBIT summary output (CSV)")
    parser.add_argument("--subtype", help="subtyper.py summary output (CSV)")
    parser.add_argument("--denovo", help="De novo assembly (FASTA)")
    parser.add_argument("--min_ncbi_stats_n", type=int, default=10,
                        help="Minimum N in NCBI stats for z-scores")
    parser.add_argument("--min_depth", type=int, default=20,
                        help="Minimum read depth for auto QC")
    parser.add_argument("--min_qual", type=float, default=0.8,
                        help="Minimum Q30 rate after filtering for auto QC")
    parser.add_argument("--max_z_score", type=float, default=2.58,
                        help="Maximum absolute z-score for auto QC")
    parser.add_argument("--log-level",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help="Logging level")
    parser.add_argument("--version", action="version",
                        version=version, help="Show script version and exit.")
    parser.add_argument("--log-file", help="Log file path (optional)")

    args = parser.parse_args()

    # configure logging
    logging.basicConfig(
        filename=args.log_file,
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    data: Dict[str, Any] = {'sample': args.sample}

    # Samplesheet (manual overrides)
    samplesheet_data = load_csv(args.samplesheet, 'samplesheet')
    data |= parse_samplesheet(samplesheet_data, args.sample)

    # Species
    species_cols = ['species', 'species_confidence']
    taxid: Optional[str] = None
    if not all(c in data for c in species_cols):
        if args.species:
            species_data = load_csv(args.species, 'species')
            if species_data:
                species_parsed, taxid = parse_species(species_data[0])
                data |= species_parsed
                species_cols = list(species_parsed.keys())
            else:
                logging.warning("Species CSV is empty")
        else:
            logging.warning("No species data provided")

    # Subtype
    subtype_cols = ['subtype']
    if not all(c in data for c in subtype_cols):
        if args.subtype:
            subtype_data = load_csv(args.subtype, 'subtype')
            if subtype_data:
                subtype_parsed = parse_subtype(subtype_data[0])
                data |= subtype_parsed
                subtype_cols = list(subtype_parsed.keys())
            else:
                logging.warning("Subtype CSV is empty")
        else:
            logging.warning("No subtype data provided")

    # De novo / assembly stats
    denovo_cols: List[str] = []
    if args.denovo:
        with load_fasta(args.denovo) as fasta_records:
            genome_stats_data = genome_stats(fasta_records)
        data |= genome_stats_data
        denovo_cols = list(genome_stats_data.keys())
    else:
        logging.warning("No assembly data provided")

    # Read stats
    read_cols: List[str] = []
    if args.read_stats:
        stats = load_json(args.read_stats, 'read_stats')
        stats_parsed = parse_read_stats(stats)
        data |= stats_parsed
        read_cols = list(stats_parsed.keys())
    else:
        logging.warning("No read stats data provided")

    # Depth estimate
    depth_col: List[str] = []
    bases = data.get('total_bases_after_filtering', 0)  # from fastp
    length = data.get('denovo_length', 0)
    if bases and length:
        try:
            b = float(bases)
            L = float(length)
            data['denovo_depth'] = int(round(b / L)) if L > 0 else 0
            depth_col = ['denovo_depth']
        except (TypeError, ValueError):
            logging.warning("Could not compute read depth due to invalid values")
    else:
        logging.warning("Read depth not estimated")

    # Taxid z-scores
    taxid_cols: List[str] = []
    if all(c in data for c in ['denovo_length', 'denovo_gc']) and args.ncbi_stats:
        ncbi_stats_data = load_csv(args.ncbi_stats, 'ncbi_stats')
        taxid_stats_data = taxid_stats(data, ncbi_stats_data, taxid, args.min_ncbi_stats_n)
        data |= taxid_stats_data
        taxid_cols = list(taxid_stats_data.keys())
    else:
        logging.warning("Z-scores not calculated for assembly length and GC")

    # Auto QC (z-scores are optional; included if present)
    data = auto_qc(data, args.min_depth, args.min_qual, args.max_z_score)

    # Normalize values
    for k, v in list(data.items()):
        # list values to semicolon-joined strings (e.g., qc_reason)
        if isinstance(v, list):
            data[k] = ';'.join(map(str, v))
        # round floats but keep depth as int
        elif isinstance(v, float) and not k.endswith('_depth'):
            data[k] = round(v, 2)

    # Output
    header_order = (['sample', 'qc_status', 'qc_reason']
                    + species_cols + subtype_cols
                    + denovo_cols + depth_col + taxid_cols + read_cols)
    data = {col: data.get(col) for col in header_order}
    with open(f"{args.sample}-summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)


if __name__ == '__main__':
    main()
