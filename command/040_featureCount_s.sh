#!/bin/sh
# featureCount for single end

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check
# output directory
OUTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/${TYPE}/${TARGET}"

# gtf file
GTFFILE="/Users/Emma/Documents/Bioinformatics/RefGenome/ENSEMBLE/Homo_sapiens.GRCh38.104.gtf"

# bam directory
BAMDIR="/Volumes/HDD24TB/${TARGET}/fastp/bam"
BAMFILES=`ls ${BAMDIR}`

cd ${BAMDIR}
# -s: # 0 (unstranded, default), 1 (stranded) and 2 (reversely stranded)
featureCounts -T 4 -t exon -g gene_id -O -M -s 2 \
-a ${GTFFILE} -o ${OUTDIR}/counts.fastp.txt ${BAMFILES}
echo finished
