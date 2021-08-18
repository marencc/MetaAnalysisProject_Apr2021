#!/bin/bash
# creating fastQC repots

# targets
TYPE="aging"  # check
TARGET="GSE129643"  # check
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}

# fastq directory to be FastQC  ## check
FILEDIR=${TARGETDIR}/fastq
# FILEDIR=${TARGETDIR}/fastp
# FILEDIR=${TARGETDIR}/sickle

# fastqc directory
mkdir -p ${TARGETDIR}/fastqc
OUTDIR=${TARGETDIR}/fastqc

# file ids
FILEIDS=${PROJECTDIR}/${TYPE}/${TARGET}/qclist.txt  # SRR_Acc_List.txt

# run FastQC
cd ${FILEDIR}
cat $FILEIDS | while read line; do
    echo FastQC: ${line}
    fastqc -t 8 --nogroup -o ${OUTDIR} ${line}*.fastq.gz
    echo number of FastQC files: `ls ${OUTDIR}/*fastqc.zip | wc -l`
    echo "\n" 
done

echo finished


# FastQC for one file
# fastqc -t 10 -o /Volumes/HDD24TB/GSE113165/fastqc /Volumes/HDD24TB/GSE113165/SRR7007949.fastq.gz