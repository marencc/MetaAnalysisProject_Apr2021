#!/bin/sh
# featureCount for single end

# targets
TYPE="aging"  # check
TARGET="GSE129643"  # check
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}

# gtf file
GTFFILE=/Volumes/HDD24TB/RefGenome/ENSEMBLE/Homo_sapiens.GRCh38.104.gtf

# bam directory
BAMDIR=${TARGETDIR}/star
BAMFILES=`ls ${BAMDIR}/*.bam`

cd ${BAMDIR}
# -s: # 0 (unstranded, default), 1 (stranded) and 2 (reversely stranded)
featureCounts -T 10 -t exon -g gene_id -O -M --fraction \
-a ${GTFFILE} -o ${TARGETDIR}/counts.txt ${BAMFILES}
echo finished
