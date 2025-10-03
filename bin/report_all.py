#!/usr/bin/env python3

import csv
import sys
import os
from pathlib import Path
from collections import OrderedDict, Counter, defaultdict
from statistics import median

def _sniff_delimiter(csvfile):
    # Robust delimiter sniff with fallback
    sample = csvfile.read(2048)
    csvfile.seek(0)
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(sample)
        return dialect.delimiter
    except Exception:
        return ','

def _consensus_field_order(headers_list):
    """
    Build a consensus column order:
      - Start with the most frequent full header sequence (mode)
      - For any remaining columns, order by median position across files,
        tie-break by frequency (desc), then column name.
    """
    if not headers_list:
        return []

    # 1) Most-common full header sequence
    header_counter = Counter(tuple(h) for h in headers_list if h)
    mode_header, _ = header_counter.most_common(1)[0]
    base_order = list(mode_header)

    # 2) Stats for median position & presence frequency
    pos_lists = defaultdict(list)   # col -> [positions...]
    freq = Counter()                # col -> in how many files it appears

    for hdr in headers_list:
        idx_map = {c: i for i, c in enumerate(hdr)}
        for c in hdr:
            pos_lists[c].append(idx_map[c])
            freq[c] += 1

    # All columns seen
    all_cols = set().union(*[set(h) for h in headers_list])

    # Remaining columns (not already in base)
    remaining = [c for c in all_cols if c not in base_order]

    # Sort remaining by median position, then by frequency desc, then name
    def sort_key(c):
        med = median(pos_lists[c]) if pos_lists[c] else float('inf')
        return (med, -freq[c], c)

    remaining_sorted = sorted(remaining, key=sort_key)

    return base_order + remaining_sorted

def combine_csv_files(*file_paths, output_file="combined_data.csv", add_source_file=False):
    """
    Combine multiple CSV files based on common column names, preserving
    the most common/consensus column order across inputs.
    """
    if not file_paths:
        print("No CSV files provided!")
        return None

    headers_list = []
    existing_files = []

    # First pass: collect headers from each file
    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"Warning: File '{file_path}' not found. Skipping...")
            continue
        try:
            with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
                delimiter = _sniff_delimiter(csvfile)
                reader = csv.DictReader(csvfile, delimiter=delimiter)
                if not reader.fieldnames:
                    print(f"Warning: No header found in '{file_path}'. Skipping...")
                    continue
                headers_list.append(reader.fieldnames)
                existing_files.append(file_path)
        except Exception as e:
            print(f"Error reading {file_path} for header detection: {e}")

    if not headers_list:
        print("No readable CSV headers found!")
        return None

    # Build consensus column order
    consensus = _consensus_field_order(headers_list)

    # Optionally add a source_file column (at the end by default)
    if add_source_file and 'source_file' not in consensus:
        consensus = consensus + ['source_file']

    print(f"Consensus column count: {len(consensus)}")
    print(f"Consensus order: {consensus}")

    # Second pass: read and combine in that order
    all_data = []
    for file_path in existing_files:
        try:
            with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
                delimiter = _sniff_delimiter(csvfile)
                reader = csv.DictReader(csvfile, delimiter=delimiter)
                file_data = []
                for row in reader:
                    new_row = OrderedDict()
                    for col in consensus:
                        if col == 'source_file':
                            new_row[col] = Path(file_path).name
                        else:
                            new_row[col] = row.get(col, '')
                    file_data.append(new_row)
                all_data.extend(file_data)
                print(f"Loaded {file_path}: {len(file_data)} rows")
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    if not all_data:
        print("No valid data was loaded!")
        return None

    # Write output
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=consensus, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_data)
        print(f"\nCombined dataset: {len(all_data)} rows, {len(consensus)} columns")
        print(f"Combined data saved to: {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        return None

    return all_data, consensus

def analyze_missing_data(data, columns):
    if not data:
        return
    print("\nMissing data analysis:")
    missing_counts = {col: 0 for col in columns}
    for row in data:
        for col in columns:
            if not (row.get(col, '') or '').strip():
                missing_counts[col] += 1
    total_rows = len(data)
    for col, count in missing_counts.items():
        if count > 0:
            percentage = (count / total_rows) * 100
            print(f"  {col}: {count}/{total_rows} ({percentage:.1f}%) missing")

def preview_data(data, num_rows=1):
    if not data:
        return
    print(f"\nData preview {min(num_rows, len(data))} rows:")
    for i, row in enumerate(data[:num_rows]):
        print(f"Row {i+1}:")
        for key, value in row.items():
            print(f"  {key}: {value}")
        print()

def main():
    if len(sys.argv) < 2:
        print("Usage: python combine_csv.py file1.csv file2.csv [file3.csv ...]")
        print("       python combine_csv.py file1.csv file2.csv --output combined.csv")
        print("Flags: --source-file     (append a source_file column)")
        return

    args = sys.argv[1:]
    output_filename = "CorgiSNPs-summary.csv"
    add_source_file = False

    # Parse flags
    if "--output" in args:
        i = args.index("--output")
        if i + 1 < len(args):
            output_filename = args[i + 1]
            args = args[:i] + args[i + 2:]

    if "--source-file" in args:
        add_source_file = True
        args = [a for a in args if a != "--source-file"]

    csv_files = args

    combined = combine_csv_files(*csv_files, output_file=output_filename, add_source_file=add_source_file)
    if combined:
        combined_data, columns = combined
        preview_data(combined_data)
        analyze_missing_data(combined_data, columns)

if __name__ == "__main__":
    main()
