#!/usr/bin/env bash
set -euo pipefail

ref_txome=$1
fastq=$2
out_bam=$3

mkdir -p "$(dirname "$out_bam")"

#minimap2 -ax map-ont "$ref_txome" "$fastq" | \
#awk '{ if($1 ~ /^@/){print $0} else {match($1, /([ACGT]+)_([ACGT]+)#/, a); if(length(a[1])>0 && length(a[2])>0) $0=$0"\tCB:Z:"a[1]"\tUB:Z:"a[2]; print $0} }' | \
#samtools sort -n -o "$out_bam"
minimap2 -ax map-ont "$ref_txome" "$fastq" | \
awk '{ if($1 ~ /^@/){print $0} else {split($1, p, "_"); if(length(p)>1){split(p[2], q, "#"); $0=$0"\tCB:Z:"p[1]"\tUB:Z:"q[1]} print $0} }' | \
samtools sort -n -o "$out_bam"

