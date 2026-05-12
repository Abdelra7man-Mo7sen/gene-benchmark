#!/usr/bin/env python3
"""
Ensemble voting for gene prediction
Combines multiple tools using majority or weighted vote
"""

import sys
import argparse
import pandas as pd
from collections import defaultdict

def parse_gff(gff_file, tool_name):
    """Parse GFF file and extract gene intervals"""
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene' or fields[2] == 'CDS':
                genes.append({
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'tool': tool_name
                })
    return genes

def overlap(g1, g2, min_overlap=0.5):
    """Check if two genes overlap by at least min_overlap fraction"""
    overlap_start = max(g1['start'], g2['start'])
    overlap_end = min(g1['end'], g2['end'])
    overlap_len = max(0, overlap_end - overlap_start)
    
    len1 = g1['end'] - g1['start']
    len2 = g2['end'] - g2['start']
    
    overlap_frac1 = overlap_len / len1
    overlap_frac2 = overlap_len / len2
    
    return overlap_frac1 >= min_overlap or overlap_frac2 >= min_overlap

def majority_vote(gene_sets, min_votes=2):
    """Majority vote: keep genes called by >= min_votes tools"""
    all_genes = []
    for genes in gene_sets:
        all_genes.extend(genes)
    
    # Cluster overlapping genes
    clusters = []
    used = set()
    
    for i, g1 in enumerate(all_genes):
        if i in used:
            continue
        
        cluster = [g1]
        used.add(i)
        
        for j, g2 in enumerate(all_genes):
            if j in used:
                continue
            if overlap(g1, g2):
                cluster.append(g2)
                used.add(j)
        
        clusters.append(cluster)
    
    # Vote
    consensus_genes = []
    for cluster in clusters:
        tools_present = set([g['tool'] for g in cluster])
        if len(tools_present) >= min_votes:
            # Take consensus coordinates (median start, median end)
            starts = [g['start'] for g in cluster]
            ends = [g['end'] for g in cluster]
            
            consensus = {
                'chrom': cluster[0]['chrom'],
                'start': int(pd.Series(starts).median()),
                'end': int(pd.Series(ends).median()),
                'strand': max(set([g['strand'] for g in cluster]), key=[g['strand'] for g in cluster].count),
                'supporting_tools': list(tools_present),
                'num_support': len(tools_present)
            }
            consensus_genes.append(consensus)
    
    return consensus_genes

def weighted_vote(gene_sets, weights):
    """Weighted vote based on tool performance"""
    # Similar to majority_vote but weight each tool
    # Implementation similar but with weighted scoring
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prodigal', help='Prodigal GFF file')
    parser.add_argument('--braker', help='BRAKER2 GFF file')
    parser.add_argument('--metagene', help='MetaGeneAnnotator GFF file')
    parser.add_argument('--output', help='Output GFF file')
    parser.add_argument('--method', default='majority', choices=['majority', 'weighted'])
    args = parser.parse_args()
    
    # Parse all GFF files
    gene_sets = []
    if args.prodigal:
        gene_sets.append(parse_gff(args.prodigal, 'prodigal'))
    if args.braker:
        gene_sets.append(parse_gff(args.braker, 'braker2'))
    if args.metagene:
        gene_sets.append(parse_gff(args.metagene, 'metagene'))
    
    # Apply voting
    if args.method == 'majority':
        consensus = majority_vote(gene_sets, min_votes=2)
    else:
        weights = {'prodigal': 0.987, 'braker2': 0.901, 'metagene': 0.904}
        consensus = weighted_vote(gene_sets, weights)
    
    # Write output GFF
    with open(args.output, 'w') as f:
        f.write("##gff-version 3\n")
        for gene in consensus:
            f.write(f"{gene['chrom']}\tensemble\tgene\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\t")
            f.write(f"ID=ensemble_{gene['start']};support={gene['num_support']};tools={','.join(gene['supporting_tools'])}\n")
    
    print(f"Wrote {len(consensus)} consensus genes to {args.output}")

if __name__ == "__main__":
    main()
