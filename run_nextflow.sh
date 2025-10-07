#!/bin/bash

module load nextflow/25.04.6 
module load apptainer

# Tell Nextflow how many SLURM CPUs
export SLURM_CPUS_PER_TASK=4

# Set default CPUs for processes
export NXF_DEFAULT_CPUS=4

# Limit Nextflow Java memory usage
export NXF_OPTS='-Xms1g -Xmx4g'


DATA_DIR=/project/6007998/maposto/PROJECTS/MPAQT_REBUTTAL/LUAD_SCLC-A_bulk_sc/DATA/ONT/SRR30947479

#  --ref_genome 'ref/genome_chr1.fa' \
#  --container '/project/6007998/_shared/apptainer/lr_sc_pipeline/lr_sc_pipeline.sif' \
#nextflow run main.nf -with-apptainer -resume \
nextflow run main.nf \
  -with-apptainer /project/6007998/_shared/apptainer/lr_sc_pipeline/lr_sc_pipeline.sif \
  -resume \
  --input "${DATA_DIR}/SRR30947479.10000_subset.fastq" \
  --outdir 'RUNS/test3' \
  --ref_txome 'ref/gencode.v38.transcripts.first1000.fa'
