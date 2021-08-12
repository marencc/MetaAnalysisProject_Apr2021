#!/bin/bash
# creating fastQC repots

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check

# fastq directory
FASTQDIR="/Volumes/HDD24TB/${TARGET}"

# sickle directory
SICKLEDIR="trimmed"
# fastp directory
FASTPDIR="fastp"


# fastqc directory
mkdir -p ${FASTQDIR}/fastqc
OUTDIR="fastqc"

# file ids
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021"
FILEIDS="${PROJECTDIR}/${TYPE}/${TARGET}/quality_check_list.txt" # SRR_Acc_List.txt


# run FastQC
cd ${FASTQDIR}
cat $FILEIDS | while read line; do
    echo FastQC: ${line}
    # check which type of files to use
    # sickle
    # fastqc -t 4 --nogroup -o ${OUTDIR} ${SICKLEDIR}/${line}.trimmed.fastq.gz
    
    # fastp
    fastqc -t 4 --nogroup -o ${OUTDIR} ${FASTPDIR}/${line}.fastp.fastq.gz
    echo number of FastQC files: `ls ${OUTDIR}/*fastp_fastqc.zip | wc -l`
    echo "\n" 
done

echo finished


# fastqc -t 4 -o /Volumes/HDD24TB/GSE113165/fastqc /Volumes/HDD24TB/GSE113165/SRR7007949.fastq.gz