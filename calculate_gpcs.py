#!/usr/bin/env python3
"""
Calculate Gene Prediction Confidence Score (GPCS)
Based on RNA-seq, proteomics, synteny, and start codon context
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import pysam

def calculate_gpcs(gene, genome_fasta, rna_bam, proteomics_file, synteny_file):
    """
    GPCS = 0.4*RNA + 0.3*Prot + 0.2*Synt + 0.1*Start
    
    Args:
        gene: dict with keys: 'chrom', 'start', 'end', 'strand'
        genome_fasta: path to genome FASTA
        rna_bam: path to RNA-seq BAM file
        proteomics_file: path to peptide evidence
        synteny_file: path to synteny conservation scores
    """
    
    # 1. RNA-seq coverage score
    rna_score = calc_rna_coverage(gene, rna_bam)
    
    # 2. Proteomic support
    prot_score = calc_proteomic_support(gene, proteomics_file)
    
    # 3. Synteny conservation
    synt_score = calc_synteny_conservation(gene, synteny_file)
    
    # 4. Start codon context
    start_score = calc_start_context(gene, genome_fasta)
    
    # Weighted sum
    gpcs = (0.4 * rna_score) + (0.3 * prot_score) + (0.2 * synt_score) + (0.1 * start_score)
    
    return gpcs

def calc_rna_coverage(gene, bam_file):
    """Calculate mean RNA-seq RPKM over gene region"""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Get coverage
    coverage_array = np.zeros(gene['end'] - gene['start'] + 1)
    for pileupcolumn in samfile.pileup(gene['chrom'], gene['start'], gene['end']):
        coverage_array[pileupcolumn.pos - gene['start']] = pileupcolumn.nsegments
    
    mean_cov = np.mean(coverage_array)
    rpkm = (mean_cov * 1000) / (len(coverage_array) / 1000)
    
    # Score: 1 if RPKM > 5, else RPKM/5
    score = min(1.0, rpkm / 5.0)
    return score

def calc_proteomic_support(gene, peptide_file):
    """Check if unique peptides map to this gene"""
    peptides = pd.read_csv(peptide_file, sep='\t')
    
    # Filter peptides mapping to this gene
    gene_peptides = peptides[
        (peptides['chrom'] == gene['chrom']) &
        (peptides['start'] >= gene['start']) &
        (peptides['end'] <= gene['end'])
    ]
    
    unique_peptides = gene_peptides[gene_peptides['unique'] == True]
    
    if len(unique_peptides) >= 2:
        return 1.0
    elif len(unique_peptides) == 1:
        return 0.5
    else:
        return 0.0

def calc_synteny_conservation(gene, synteny_file):
    """Fraction of gene aligned in orthologous genomes"""
    synteny_data = pd.read_csv(synteny_file, sep='\t')
    gene_synt = synteny_data[
        (synteny_data['gene_id'] == gene['id'])
    ]
    
    if len(gene_synt) == 0:
        return 0.0
    
    conservation = gene_synt['aligned_fraction'].mean()
    return min(1.0, conservation)

def calc_start_context(gene, genome_fasta):
    """Evaluate ribosome binding site (prokaryote) or Kozak (eukaryote)"""
    # Extract sequence around start codon (-20 to +10)
    record = SeqIO.read(genome_fasta, "fasta")
    
    if gene['strand'] == '+':
        start = gene['start'] - 21
        end = gene['start'] + 10
        seq = str(record.seq[start:end])
    else:
        # Reverse complement for minus strand
        start = gene['end'] - 10
        end = gene['end'] + 21
        seq = str(record.seq[start:end].reverse_complement())
    
    # Check for Shine-Dalgarno (AGGAGG) in prokaryotes
    if 'AGGAGG' in seq[-20:-5]:
        return 1.0
    elif 'GGAGG' in seq[-20:-5]:
        return 0.75
    elif 'GGA' in seq[-20:-5]:
        return 0.5
    else:
        return 0.25

def classify_gpcs(gpcs):
    """Classify confidence level"""
    if gpcs >= 0.6:
        return "HIGH"
    elif gpcs >= 0.3:
        return "MODERATE"
    else:
        return "LOW"

if __name__ == "__main__":
    import os
    
    gene_example = {
        'id': 'gene_001',
        'chrom': 'NC_000913.3',
        'start': 1000,
        'end': 2000,
        'strand': '+'
    }
    
    # Check if files exist before running
    fasta_file = "data/E_coli.fasta"
    bam_file = "data/rna_seq.bam"
    peptides_file = "data/peptides.tsv"
    synteny_file = "data/synteny.tsv"
    
    if not os.path.exists(fasta_file):
        print(f"Warning: {fasta_file} not found. Using dummy data.")
        gpcs = 0.5
    else:
        gpcs = calculate_gpcs(
            gene_example,
            fasta_file,
            bam_file if os.path.exists(bam_file) else None,
            peptides_file if os.path.exists(peptides_file) else None,
            synteny_file if os.path.exists(synteny_file) else None
        )
    
    print(f"GPCS: {gpcs:.3f} - {classify_gpcs(gpcs)}")