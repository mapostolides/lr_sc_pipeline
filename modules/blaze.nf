process blaze {
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq
    path known_bc_list

    output:
    path "blaze/matched_reads.fastq.gz", emit: filtered_fastq
    path "blaze/*"

    script:
    def bc_arg = known_bc_list.name != 'NO_FILE' ? "--known-bc-list ${known_bc_list}" : ""
    def msg = known_bc_list.name != 'NO_FILE' ? "Using whitelist: ${known_bc_list}" : "No whitelist provided"
    """
    echo "[INFO] ${msg}"   
    blaze --expect-cells 5000 ${bc_arg} \
          --output-prefix blaze/ \
          --threads 4 \
          $fastq
    """
}