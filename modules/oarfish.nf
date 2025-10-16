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