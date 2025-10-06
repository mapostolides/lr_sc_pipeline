nextflow.enable.dsl=2
params.input  = "SRR30947479.10000_subset.fastq"
params.outdir = "RUNS/test2"
params.ref_genome    = "ref/genome_chr1.fa"
params.ref_txome    = "ref/gencode.v38.transcripts.first1000.fa"

process blaze {
    container 'lr_sc_pipeline:dev'
    
	//  put actual outfiles inside specified directory
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq

    output:
	path "blaze/matched_reads.fastq.gz", emit: filtered_fastq

    script:
    """
    blaze --expect-cells 5000 \
          --output-prefix blaze/ \
          --threads 1 \
          $fastq
    """
}

process minimap2_genome {
    container 'lr_sc_pipeline:dev'
    publishDir params.outdir, mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    path fastq
	path ref_genome

    output:
    path "minimap2/*"

    script:
    """
    mkdir -p minimap2
    minimap2 -ax splice $ref_genome $fastq | samtools view -bS - > minimap2/aligned.bam 
    """
}

process minimap2_txome {
    container 'lr_sc_pipeline:dev'
    publishDir params.outdir, mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    path fastq
    path ref_txome

    output:
    path "minimap2/*"

    script:
    """
    mkdir -p minimap2
	minimap2 -ax map-ont $ref_txome $fastq | samtools view -bS - > minimap2/aln.txome.bam

    """
}

workflow {
    fastq = file(params.input)
	ref_genome = file(params.ref_genome)
	ref_txome = file(params.ref_genome)
    blaze_out = blaze(fastq)
	//minimap2(blaze_out.filtered_fastq, ref_genome)
	minimap2_txome(blaze_out.filtered_fastq, ref_txome)
}
