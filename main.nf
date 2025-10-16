nextflow.enable.dsl=2

//params.container = '/project/6007998/_shared/apptainer/lr_sc_pipeline/lr_sc_pipeline.sif'
//params.input  = "SRR30947479.10000_subset.fastq"
//params.input  = "SRR30947479.1000_subset.fastq"
//params.outdir = "RUNS/test3"
//params.ref_genome    = "ref/genome_chr1.fa"
//params.ref_txome    = "ref/gencode.v38.transcripts.first1000.fa"

include { blaze } from './modules/blaze'
include { minimap2_txome } from './modules/minimap2'
include { oarfish } from './modules/oarfish'

	//ref_genome = file(params.ref_genome)
workflow {
    fastq_ch = Channel.fromPath(params.input_fastq)
	// Define optional file parameter
    if (params.known_bc_list) {
        log.info "[INFO] Using whitelist: ${params.known_bc_list}"
    } else {
        log.info "[INFO] No whitelist provided"
    }
    whitelist_ch = params.known_bc_list ? 
        Channel.fromPath(params.known_bc_list) : Channel.value(file('NO_FILE'))
    // reference transcriptome
	ref_txome_ch = Channel.fromPath(params.ref_txome)
	// STEP1: BLAZE
    blaze_out = blaze(fastq_ch, whitelist_ch)  
	// STEP2: MINIMAP
	//minimap2_out = minimap2_txome(blaze_out.filtered_fastq, ref_txome_ch)
	// STEP3: OARFISH
	//oarfish(minimap2_out.txome_bam)
}
