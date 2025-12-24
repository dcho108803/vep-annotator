#!/bin/bash
#
# VEP Annotator Data Download Script
#
# Downloads required annotation data files for full VEP functionality.
# Run this script to set up all data files.
#
# Usage: ./download_data.sh [output_dir]
#

set -e

OUTPUT_DIR="${1:-./data}"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "==================================="
echo "VEP Annotator Data Download Script"
echo "==================================="
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to download with progress
download() {
    local url="$1"
    local output="$2"
    echo "Downloading: $output"
    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "$output" "$url"
    elif command -v curl &> /dev/null; then
        curl -L -o "$output" "$url"
    else
        echo "Error: Neither wget nor curl found"
        exit 1
    fi
}

# ============================================================================
# Required: Gene annotations (GTF) and Reference genome (FASTA)
# ============================================================================

echo ""
echo "1. Downloading Ensembl GTF annotation..."
if [ ! -f "Homo_sapiens.GRCh38.110.gtf.gz" ]; then
    download "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz" \
        "Homo_sapiens.GRCh38.110.gtf.gz"
fi

echo ""
echo "2. Downloading Reference genome (primary assembly)..."
echo "   NOTE: This file is ~900MB compressed, ~3GB uncompressed"
if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then
    download "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fi

# ============================================================================
# Conservation Scores (bigWig format)
# ============================================================================

echo ""
echo "3. Downloading conservation scores..."

mkdir -p conservation
cd conservation

if [ ! -f "hg38.phyloP100way.bw" ]; then
    echo "   PhyloP 100-way (~10GB)..."
    download "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw" \
        "hg38.phyloP100way.bw"
fi

if [ ! -f "hg38.phastCons100way.bw" ]; then
    echo "   PhastCons 100-way (~10GB)..."
    download "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw" \
        "hg38.phastCons100way.bw"
fi

cd ..

# ============================================================================
# Regulatory Build (GFF3 format)
# ============================================================================

echo ""
echo "4. Downloading Ensembl Regulatory Build..."

mkdir -p regulatory
cd regulatory

if [ ! -f "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.gff.gz" ]; then
    download "https://ftp.ensembl.org/pub/release-110/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz" \
        "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.gff.gz"
fi

cd ..

# ============================================================================
# ClinVar (VCF format)
# ============================================================================

echo ""
echo "5. Downloading ClinVar..."

mkdir -p clinvar
cd clinvar

if [ ! -f "clinvar.vcf.gz" ]; then
    download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
        "clinvar.vcf.gz"
    download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" \
        "clinvar.vcf.gz.tbi"
fi

cd ..

# ============================================================================
# Print Summary
# ============================================================================

echo ""
echo "==================================="
echo "Download Summary"
echo "==================================="
echo ""
echo "Downloaded files are in: $OUTPUT_DIR"
echo ""
echo "Core files (required):"
echo "  - GTF: Homo_sapiens.GRCh38.110.gtf.gz"
echo "  - FASTA: Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
echo ""
echo "Conservation scores:"
echo "  - conservation/hg38.phyloP100way.bw"
echo "  - conservation/hg38.phastCons100way.bw"
echo ""
echo "Regulatory:"
echo "  - regulatory/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.gff.gz"
echo ""
echo "Clinical:"
echo "  - clinvar/clinvar.vcf.gz"
echo ""
echo "==================================="
echo "Manual Downloads Required"
echo "==================================="
echo ""
echo "The following files require manual download due to licensing:"
echo ""
echo "1. dbNSFP (pathogenicity scores):"
echo "   https://sites.google.com/site/jpaboratory/dbNSFP"
echo "   Download: dbNSFP4.4a.zip"
echo "   Extract and index: tabix -s 1 -b 2 -e 2 -c '#' dbNSFP4.4a.txt.gz"
echo ""
echo "2. SpliceAI (splice predictions):"
echo "   https://basespace.illumina.com/s/otSPW8hnhaZR"
echo "   Requires Illumina BaseSpace account"
echo ""
echo "3. gnomAD (population frequencies):"
echo "   https://gnomad.broadinstitute.org/downloads"
echo "   Download exomes and/or genomes VCF"
echo ""
echo "==================================="
echo "Example usage:"
echo "==================================="
echo ""
echo "vep_annotator --gtf $OUTPUT_DIR/Homo_sapiens.GRCh38.110.gtf.gz \\"
echo "              --fasta $OUTPUT_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \\"
echo "              --phylop $OUTPUT_DIR/conservation/hg38.phyloP100way.bw \\"
echo "              --phastcons $OUTPUT_DIR/conservation/hg38.phastCons100way.bw \\"
echo "              --regulatory $OUTPUT_DIR/regulatory/*.gff.gz \\"
echo "              --annotation clinvar:$OUTPUT_DIR/clinvar/clinvar.vcf.gz \\"
echo "              -v 7:140753336:A:T"
echo ""
