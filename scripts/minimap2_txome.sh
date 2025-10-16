#!/usr/bin/env bash
set -euo pipefail

ref_txome=$1
fastq=$2
out_bam=$3

mkdir -p "$(dirname "$out_bam")"



# minimap2 -ax map-ont "$ref_txome" "$fastq" | \
# awk '{ if($1 ~ /^@/){print $0} else {split($1, p, "_"); if(length(p)>1){split(p[2], q, "#"); $0=$0"\tCB:Z:"p[1]"\tUB:Z:"q[1]} print $0} }' | \
# samtools sort -n -o "$out_bam"

# --- PIPELINE SUMMARY ---
# 1. ALIGNMENT: Uses Minimap2 (map-ont) to align long reads from $fastq to the reference $ref_txome.
# 2. TAGGING (AWK): Parses the read name (expected format: BARCODE_UMI#...) to extract and inject 
#    standard single-cell tags (CB:Z: for Cell Barcode and UB:Z: for UMI) into the SAM stream.
# 3. SORTING: Uses Samtools sort -n to convert the stream to BAM and sort alignments strictly by read name, 
#    which is MANDATORY for subsequent single-cell quantification tools.
# --------------------------

minimap2 -ax map-ont "$ref_txome" "$fastq" | \
awk '{
    # Check if line is a SAM header (starts with @)
    if ($1 ~ /^@/) {
        print $0
    } 
    # Process alignment records
    else {
        # Split read name by the first underscore
        split($1, p, "_"); 
        
        # Check if split was successful (i.e., barcode and UMI parts exist)
        if (length(p) > 1) {
            # Split the second part (p[2]) to isolate the UMI from its identifier (#)
            split(p[2], q, "#"); 
            
            # Inject Cell Barcode (CB:Z) and UMI (UB:Z) tags into the SAM record ($0)
            $0 = $0 "\tCB:Z:" p[1] "\tUB:Z:" q[1]
        } 
        print $0
    } 
}' | \
samtools sort -n -o "$out_bam"




