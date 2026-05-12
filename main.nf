#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.genome_fasta = "data/input/genome.fasta"
params.genome_type = "prokaryote"  // or "eukaryote"
params.outdir = "results"
params.rna_seq_bam = false
params.threads = 8

// Define processes
process run_prodigal {
    tag "Prodigal on ${genome_id}"
    publishDir "${params.outdir}/prodigal", mode: 'copy'
    
    input:
    path genome
    
    output:
    path "prodigal.gff"
    
    script:
    """
    prodigal -i ${genome} -o prodigal.gff -f gff -p single
    """
}

process run_braker2 {
    tag "BRAKER2 on ${genome_id}"
    publishDir "${params.outdir}/braker2", mode: 'copy'
    
    input:
    path genome
    path rna_bam
    
    output:
    path "braker.gff3"
    
    when: params.genome_type == "eukaryote" && rna_bam != false
    
    script:
    """
    braker.pl --genome=${genome} --bam=${rna_bam} \
              --threads=${params.threads} --gff3
    """
}

process run_metagene {
    tag "MetaGeneAnnotator on ${genome_id}"
    publishDir "${params.outdir}/metagene", mode: 'copy'
    
    input:
    path genome
    
    output:
    path "metagene.gff"
    
    script:
    """
    metagene_annotator ${genome} > metagene.gff
    """
}

// Workflow
workflow {
    // Read genome
    genome_ch = Channel.fromPath(params.genome_fasta)
    
    // Run prokaryotic tools
    if (params.genome_type == "prokaryote") {
        prodigal_ch = run_prodigal(genome_ch)
        metagene_ch = run_metagene(genome_ch)
    }
    
    // Run eukaryotic tools
    if (params.genome_type == "eukaryote" && params.rna_seq_bam) {
        braker_ch = run_braker2(genome_ch, file(params.rna_seq_bam))
    }
    
    // Run ensemble
    ensemble_ch = run_ensemble(prodigal_ch, braker_ch, metagene_ch)
}

process run_ensemble {
    input:
    path prodigal_gff
    path braker_gff
    path metagene_gff
    
    output:
    path "ensemble.gff"
    
    when: prodigal_gff || braker_gff || metagene_gff
    
    script:
    def args = []
    if (prodigal_gff) args << "--prodigal ${prodigal_gff}"
    if (braker_gff) args << "--braker ${braker_gff}"
    if (metagene_gff) args << "--metagene ${metagene_gff}"
    
    """
    python scripts/ensemble_vote.py ${args.join(' ')} --output ensemble.gff
    """
}