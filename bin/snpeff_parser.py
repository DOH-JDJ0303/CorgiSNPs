#!/usr/bin/env python3
"""
Simplified SnpEff VCF Parser

A streamlined tool for parsing SnpEff-annotated VCF files and extracting 
variant annotations with sample genotype information.
"""

import argparse
import logging
import sys
import json
import csv
from pathlib import Path
from typing import List, Dict, Any, Optional


def setup_logging(log_level='INFO', log_file=None):
    """Setup basic logging"""
    level = getattr(logging, log_level.upper())
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    logger = logging.getLogger()
    logger.setLevel(level)
    logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)


def load_target_regions(regions_file: str) -> List[Dict]:
    """Load genomic regions from JSON file for target annotation"""
    if not regions_file:
        return []
    
    try:
        with open(regions_file, 'r') as f:
            regions_data = json.load(f)
        
        total_regions = sum(len(gene.get('regions', [])) for gene in regions_data)
        logging.info(f"Loaded {len(regions_data)} genes with {total_regions} target regions")
        return regions_data
        
    except Exception as e:
        logging.error(f"Error loading regions file: {e}")
        return []


def get_target_genes(regions_data: List[Dict]) -> List[str]:
    """Extract list of target gene names from regions data"""
    target_genes = []
    for gene_info in regions_data:
        gene_name = gene_info.get('gene', '')
        if gene_name and gene_name not in target_genes:
            target_genes.append(gene_name)
    return target_genes


def filter_variants_by_target_genes(data: List[Dict[str, Any]], target_genes: List[str]) -> List[Dict[str, Any]]:
    """Filter variants to only include those with Gene_Name matching target genes"""
    if not target_genes:
        return data
    
    filtered_data = []
    for record in data:
        gene_name = record.get('Gene_Name', '')
        if gene_name in target_genes:
            filtered_data.append(record)
    
    logging.info(f"Filtered {len(data)} records down to {len(filtered_data)} records matching target genes")
    return filtered_data


def annotate_with_targets(chrom: str, pos: int, regions_data: List[Dict]) -> Dict[str, str]:
    """Annotate variant with target gene and region information"""
    annotation = {'Target_Gene': '', 'Target_Region': ''}
    
    if not regions_data:
        return annotation
    
    try:
        pos = int(pos)
    except (ValueError, TypeError):
        logging.warning(f"Invalid position: {pos}")
        return annotation
    
    for gene_info in regions_data:
        gene_chrom = gene_info.get('chrom', '')
        gene_name = gene_info.get('gene', '')
        gene_coords = gene_info.get('coords', [])
        
        # Skip if different chromosome
        if chrom != gene_chrom:
            continue
        
        # Check if variant is within gene coordinates
        if len(gene_coords) >= 2:
            gene_start, gene_end = gene_coords[0], gene_coords[1]
            
            if gene_start <= pos <= gene_end:
                annotation['Target_Gene'] = gene_name
                
                # Check specific regions within the gene
                for region in gene_info.get('regions', []):
                    region_name = region.get('name', '')
                    region_coords = region.get('coords', [])
                    
                    if len(region_coords) >= 2:
                        region_start, region_end = region_coords[0], region_coords[1]
                        
                        if region_start <= pos <= region_end:
                            annotation['Target_Region'] = region_name
                            break
                
                # Found the gene, stop searching
                break
    
    return annotation


def parse_genotype(genotype_str: str) -> List[int]:
    """Parse genotype string to extract allele indices"""
    if not genotype_str or genotype_str in ['.', './.', '.|.']:
        return []
    
    # Handle haploid genotype
    if genotype_str.isdigit():
        return [int(genotype_str)]
    
    # Handle diploid genotype
    if '/' in genotype_str or '|' in genotype_str:
        separator = '/' if '/' in genotype_str else '|'
        alleles = genotype_str.split(separator)
        return [int(a) for a in alleles if a.isdigit()]
    
    return []


def parse_info_field(info_string: str) -> Dict[str, Any]:
    """Parse VCF INFO field into dictionary"""
    info_dict = {}
    if not info_string or info_string == '.':
        return info_dict
    
    for part in info_string.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[part] = True
    
    return info_dict


def parse_ann_field(ann_string: str, ann_fields: List[str]) -> Dict[str, List[Dict]]:
    """Parse SnpEff ANN field into structured annotations"""
    ann_dict = {}
    
    if not ann_string:
        return ann_dict
    
    annotations = ann_string.split(',')
    
    for annotation in annotations:
        parts = annotation.split('|')
        
        if len(parts) < len(ann_fields):
            logging.warning(f"Annotation has insufficient fields: {len(parts)} < {len(ann_fields)}")
            continue
        
        # Create annotation dictionary
        allele_data = {
            ann_fields[i]: parts[i] for i in range(len(ann_fields))
        }
        
        allele = allele_data.get('Allele', '')
        if allele:
            if allele not in ann_dict:
                ann_dict[allele] = []
            ann_dict[allele].append(allele_data)
    
    return ann_dict


def extract_ann_fields_from_header(header_lines: List[str]) -> List[str]:
    """Extract ANN field names from VCF header"""
    default_fields = [
        'Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID',
        'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c',
        'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length',
        'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO'
    ]
    
    for line in header_lines:
        if '##INFO=<ID=ANN' in line:
            try:
                ann_desc = line.split('Description="')[1].split('"')[0]
                
                if "Format:" in ann_desc:
                    format_part = ann_desc.split("Format:")[1].strip()
                    fields = [f.strip("' ") for f in format_part.split('|')]
                    logging.info(f"Found ANN format with {len(fields)} fields")
                    return fields
            except (IndexError, ValueError):
                logging.warning("Could not parse ANN definition from header")
    
    logging.info("Using default SnpEff annotation fields")
    return default_fields


def parse_vcf_file(filepath: str, regions_data: List[Dict]) -> List[Dict[str, Any]]:
    """Parse entire VCF file and return structured variant data"""
    if not Path(filepath).exists():
        raise FileNotFoundError(f"VCF file not found: {filepath}")
    
    logging.info(f"Parsing VCF file: {filepath}")
    
    data = []
    sample_names = []
    ann_fields = []
    header_lines = []
    variant_count = 0
    
    with open(filepath, 'r') as file:
        for line_num, line in enumerate(file, 1):
            line = line.strip()
            
            if not line:
                continue
            
            # Process header lines
            if line.startswith('#'):
                header_lines.append(line)
                
                if '#CHROM' in line:
                    # Extract column names and sample names
                    columns = line.replace('#', '').split('\t')
                    
                    try:
                        format_idx = columns.index('FORMAT')
                        sample_names = columns[format_idx + 1:]
                        logging.info(f"Found {len(sample_names)} samples")
                    except ValueError:
                        logging.warning("No FORMAT field found in VCF header")
                
                continue
            
            # Extract ANN fields from header if not already done
            if not ann_fields:
                ann_fields = extract_ann_fields_from_header(header_lines)
            
            # Parse variant line
            try:
                variant_records = parse_variant_line(
                    line, sample_names, ann_fields, regions_data
                )
                data.extend(variant_records)
                variant_count += 1
                
                if variant_count % 1000 == 0:
                    logging.info(f"Processed {variant_count} variants...")
                    
            except Exception as e:
                logging.error(f"Error parsing line {line_num}: {e}")
                continue
    
    logging.info(f"Successfully parsed {variant_count} variants, generating {len(data)} annotation records")
    return data


def parse_variant_line(line: str, sample_names: List[str], ann_fields: List[str], 
                      regions_data: List[Dict]) -> List[Dict[str, Any]]:
    """Parse a single variant line and return all relevant records"""
    fields = line.split('\t')
    
    if len(fields) < 8:  # Minimum VCF fields
        logging.warning("Skipped variant: insufficient VCF fields")
        return []
    
    # Extract basic variant information
    chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
    
    # Validate essential fields
    if not all([chrom, pos, ref, alt]):
        logging.warning(f"Skipped variant: missing essential fields")
        return []
    
    if not pos.isdigit():
        logging.warning(f"Skipped variant {chrom}: invalid position '{pos}'")
        return []
    
    # Parse INFO field
    info_dict = parse_info_field(info)
    
    # Parse ANN annotations
    ann_data = parse_ann_field(info_dict.get('ANN', ''), ann_fields)
    
    # Get target annotations
    target_annotation = annotate_with_targets(chrom, int(pos), regions_data)
    
    # Parse sample data if available
    variant_records = []
    alt_alleles = alt.split(',')
    all_alleles = [ref] + alt_alleles
    
    if len(fields) > 8 and sample_names:
        format_field = fields[8]
        format_keys = format_field.split(':')
        
        for i, sample_name in enumerate(sample_names):
            sample_field_idx = 9 + i
            if sample_field_idx >= len(fields):
                continue
            
            sample_values = fields[sample_field_idx].split(':')
            sample_info = {
                format_keys[j]: sample_values[j] if j < len(sample_values) else '.'
                for j in range(len(format_keys))
            }
            
            genotype = sample_info.get('GT', '.')
            
            # Skip samples with missing or reference-only genotypes
            if genotype in ['.', './.', '.|.']:
                continue
            
            allele_indices = parse_genotype(genotype)
            if not allele_indices or all(idx == 0 for idx in allele_indices):
                continue
            
            # Create base record for this sample
            base_record = {
                'Sample': sample_name,
                'Chrom': chrom,
                'Pos': pos,
                'Ref': ref,
                'Alt': alt,
                'Genotype': genotype,
                **{f"Sample_{key}": value for key, value in sample_info.items()},
                **target_annotation
            }
            
            # Add annotations for each non-reference allele
            sample_alleles = [all_alleles[idx] for idx in allele_indices if idx < len(all_alleles) and idx > 0]
            
            if sample_alleles:
                records_added = 0
                for allele in sample_alleles:
                    if allele in ann_data:
                        for annotation in ann_data[allele]:
                            record = base_record.copy()
                            record.update(annotation)
                            variant_records.append(record)
                            records_added += 1
                
                # If no specific annotations found, still include the base record
                if records_added == 0:
                    variant_records.append(base_record)
            
    else:
        # No sample data, just create basic variant record
        base_record = {
            'Chrom': chrom,
            'Pos': pos,
            'Ref': ref,
            'Alt': alt,
            **target_annotation
        }
        
        # Add annotations for all alt alleles
        for allele in alt_alleles:
            if allele in ann_data:
                for annotation in ann_data[allele]:
                    record = base_record.copy()
                    record.update(annotation)
                    variant_records.append(record)
            else:
                variant_records.append(base_record)
    
    return variant_records


def write_output(data: List[Dict[str, Any]], output_file: str, output_format: str) -> None:
    """Write parsed data to output file"""
    output_path = Path(output_file)
    logging.info(f"Writing {len(data)} records to {output_path} in {output_format} format")
    
    if output_format == 'json':
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    elif output_format in ['csv', 'tsv']:
        if not data:
            logging.warning("No data to write")
            return
        
        delimiter = ',' if output_format == 'csv' else '\t'
        fieldnames = data[0].keys()
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(data)
    
    logging.info(f"Output written successfully to {output_path}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Parse SnpEff-annotated VCF files and extract variant annotations"
    )
    
    parser.add_argument('input_vcf', help='Path to SnpEff-annotated VCF file')
    parser.add_argument('-o', '--outdir', default='.', help='Output directory (default: current directory)')
    parser.add_argument('-r', '--regions', help='Path to JSON file containing genomic regions')
    parser.add_argument('-f', '--format', choices=['json', 'csv', 'tsv'], 
                       default='csv', help='Output format (default: csv)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--log-file', help='Optional log file path')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level, args.log_file)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define fixed output filenames
    extension = f'.{args.format}'
    full_output = output_dir / f'variants_full{extension}'
    targets_output = output_dir / f'variants_targets{extension}'
    
    try:
        # Load target regions if provided
        regions_data = load_target_regions(args.regions)
        
        # Parse VCF file (full output)
        data = parse_vcf_file(args.input_vcf, regions_data)
        
        # Write full output
        write_output(data, str(full_output), args.format)
        
        # Create filtered output if target genes are available
        if regions_data:
            target_genes = get_target_genes(regions_data)
            filtered_data = filter_variants_by_target_genes(data, target_genes)
            write_output(filtered_data, str(targets_output), args.format)
            logging.info(f"Created: {full_output} ({len(data)} records), {targets_output} ({len(filtered_data)} records)")
        else:
            logging.info(f"Created: {full_output} ({len(data)} records)")
        
    except Exception as e:
        logging.error(f"Error during VCF parsing: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()