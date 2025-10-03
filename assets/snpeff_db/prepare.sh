#!/usr/bin/env bash
set -euo pipefail

ACCESSION="${1}"
NAME="${2}"

DATA_DIR="data/${NAME}"
mkdir -p "${DATA_DIR}"
pushd "${DATA_DIR}" >/dev/null

# Download from NCBI Datasets
datasets download genome accession "${ACCESSION}" \
  --include gff3,genome \
  --filename ncbi_dataset.zip

unzip -o ncbi_dataset.zip

cat ncbi_dataset/data/*/*.fna > sequences.fa
cat ncbi_dataset/data/*/*.gff > genes.gff

# Normalize vertical bars in BOTH GFF and FASTA so IDs stay in sync
# (SnpEff uses '|' as ANN field delimiter)
sed -i 's/|/_/g' genes.gff
# Replace only in FASTA headers (not sequences)
awk 'BEGIN{OFS=""} /^>/ {gsub(/\|/,"_"); print; next} {print}' sequences.fa > sequences.tmp && mv sequences.tmp sequences.fa

# Optional: tidy up
rm -rf ncbi_dataset* README.md md5sum.txt ncbi_dataset/data/*/README* || true

# Derive CDS/protein for your own use (SnpEff only needs genes.gff + sequences.fa)
gffread genes.gff -g sequences.fa -x cds.fa -y protein.fa

popd >/dev/null

# Update snpEff.config (add a blank line + entry)
# Format: <genome_id>.genome : <description>
# We’ll set description = NAME for simplicity
printf "\n# %s\n%s.genome : %s\n" "${NAME}" "${NAME}" "${NAME}" >> snpEff.config

# Build database from our data directory
snpEff build -c snpEff.config -dataDir data/ -gff3 -v "${NAME}"
