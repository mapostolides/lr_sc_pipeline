process minimap2_txome {
//    container params.container 
    publishDir params.outdir, mode: 'copy'
    cpus 4
    //memory '8 GB'

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