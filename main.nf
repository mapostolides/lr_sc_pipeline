nextflow.enable.dsl=2
params.input  = "SRR30947479.10000_subset.fastq"
params.outdir = "RUNS/test2"

process blaze {
    container 'lr_sc_pipeline:dev'
    
	//  put actual outfiles inside specified directory
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq

    output:
    path "blaze_out/*" 

    script:
    """
    blaze --expect-cells 5000 \
          --output-prefix blaze_out/ \
          --threads 1 \
          $fastq
    """
}

workflow {
    fastq = file(params.input)
    blaze(fastq)
}

