nextflow.enable.dsl=2
//params.container = '/project/6007998/_shared/apptainer/lr_sc_pipeline/lr_sc_pipeline.sif'
//params.input  = "SRR30947479.10000_subset.fastq"
//params.input  = "SRR30947479.1000_subset.fastq"
//params.outdir = "RUNS/test3"
//params.ref_genome    = "ref/genome_chr1.fa"
//params.ref_txome    = "ref/gencode.v38.transcripts.first1000.fa"

process blaze {
//    container params.container 
    
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
//    container params.container 
    publishDir params.outdir, mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    path fastq
	path ref_genome

    output:
	path "minimap2/aln.genome.bam", emit: genome_bam

    script:
    """
    mkdir -p minimap2
    """
    //minimap2 -ax splice $ref_genome $fastq | samtools view -bS - > minimap2/aln.genome.bam 
    //minimap2 -ax splice $ref_genome blaze/matched_reads.fastq.gz  | samtools view -bS - > minimap2/aln.genome.bam 

}

process minimap2_txome {
//    container params.container 
    publishDir params.outdir, mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    path fastq
    path ref_txome

    output:
	path "minimap2/aln.txome.name_sorted.bam", emit: txome_bam

	//minimap2 -ax map-ont $ref_txome $fastq | samtools view -bS - > minimap2/aln.txome.bam
	//samtools sort -n minimap2/aln.txome.bam > minimap2/aln.txome.sorted.bam
	//samtools index minimap2/aln.txome.sorted.bam
	//minimap2 -ax map-ont $ref_txome $fastq |  samtools sort -n -o minimap2/aln.txome.name_sorted.bam 
    script:
    """ 
    mkdir -p minimap2
	${projectDir}/scripts/minimap2_txome.sh $ref_txome $fastq minimap2/aln.txome.name_sorted.bam
    """
	//${projectDir}/scripts/minimap2_txome.sh $ref_txome blaze/matched_reads.fastq.gz minimap2/aln.txome.name_sorted.bam


}

process oarfish {
//  container params.container
  publishDir params.outdir, mode: 'copy'
  //memory '8 GB'

  input:
    path txome_bam

  output:
    path "oarfish/*"

  script:
  """
  oarfish --single-cell -j 4 -o oarfish_output -a $txome_bam
  mkdir -p oarfish; mv oarfish_output* oarfish/
  """
}

	//ref_genome = file(params.ref_genome)
	//minimap2(blaze_out.filtered_fastq, ref_genome)
workflow {
    fastq = file(params.input)
	ref_txome = file(params.ref_txome)
    blaze_out = blaze(fastq)
	minimap2_out = minimap2_txome(blaze_out.filtered_fastq, ref_txome)
	oarfish(minimap2_out.txome_bam)
}
