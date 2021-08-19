#!bin/sh
# checking strandness

# targets
TYPE="aging"  # check
TARGET="GSE129643"  # check
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}  # output directory

# input: bam file directory 
BAMDIR=${TARGETDIR}/star

# ref: bed file directory
BEDFILE=/Volumes/HDD24TB/RefGenome/ENSEMBLE/Homo_sapiens.GRCh38.104.bed

# infer_experiment.py
for id in `cat ${TARGETDIR}/qclist.txt | head -1`; do
    echo checking strandness: ${id} "(${TARGET})"
    echo "${id}.Aligned.sortedByCoord.out.bam" > ${TARGETDIR}/strandstat.txt
    
    infer_experiment.py \
    -r ${BEDFILE} \
    -i ${BAMDIR}/${id}.Aligned.sortedByCoord.out.bam \
    >> ${TARGETDIR}/strandstat.txt
done
