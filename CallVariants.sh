#!/bin/bash

# Reference genome (update this if needed)
REFERENCE="Homo_sapiens_assembly19.fixed.fasta"

# Output directory for VCF files
OUTDIR="VCF_Files"
mkdir -p "$OUTDIR"

# Loop through all BAM files in the directory
for BAM in *.bam; do

    SAMPLE_NAME=$(basename "$BAM" .bam)
    VCF_FILE="$OUTDIR/${SAMPLE_NAME}.vcf.gz"

    echo "🔄 Processing $SAMPLE_NAME..."

    # Remove any existing VCF file to prevent old data issues (optional)
    if [[ -f "$VCF_FILE" ]]; then
        echo "🧹 Removing existing VCF file: $VCF_FILE"
        rm "$VCF_FILE"
    fi

    # Re-run variant calling, even if VCF file exists
    echo "🧬 Calling variants for $SAMPLE_NAME..."
    bcftools mpileup -Ou -f "$REFERENCE" "$BAM" | \
    bcftools call -mv -Oz -o "$VCF_FILE"
    
    # Index the new VCF file
    bcftools index "$VCF_FILE"
    echo "📄 VCF created: $VCF_FILE"
done

echo "🎉 Variant calling completed for all BAM files!"
