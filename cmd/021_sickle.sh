#!/bin/sh
# trimming low quality reads by sickle

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check

# fastq directory
FASTQDIR="/Volumes/HDD24TB/${TARGET}"

# output directory
mkdir -p ${FASTQDIR}/sickle
OUTDIR="sickle"

# file ids
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021"
FILEIDS="${PROJECTDIR}/${TYPE}/${TARGET}/SRR_Acc_List.txt"

cd ${FASTQDIR} 
cat $FILEIDS | while read line; do 
    echo sickle: ${line} "(${TARGET})"
    sickle se -f ${line}.fastq.gz -t sanger -o ${OUTDIR}/${line}.sickle.fastq.gz -g 
    echo "\n" 
done 

echo finished 

