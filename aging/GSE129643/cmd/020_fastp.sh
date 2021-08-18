#!/bin/sh
# trimming low quality reads by fastp
# adapter trimming is enabled by default
# detecting adapter is TRUE for single end by default

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check

# fastq directory
FASTQDIR="/Volumes/HDD24TB/${TARGET}"

# output directory
mkdir -p ${FASTQDIR}/fastp
FASTPDIR="${FASTQDIR}/fastp"

# file ids
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021"
FILEIDS="${PROJECTDIR}/${TYPE}/${TARGET}/SRR_Acc_List.txt"

cd ${FASTPDIR}
mkdir -p html json
cat $FILEIDS | while read line; do
    echo fastp: ${line}
    # single end
    # -3: enable per read cutting by quality in tail (3'), default is disabled
    # (WARNING: this will interfere deduplication for SE data)
    fastp -i ${FASTQDIR}/${line}.fastq.gz -o ${line}.fastp.fastq.gz \
    -h html/${line}.report.html -j json/${line}.report.json --thread 4 -q 20 
    
    # pair end
    # fastp -i ${FASTQDIR}/${line}_1.fastq.gz -I ${FASTQDIR}/${line}_2.fastq.gz -3 \
    # -o ${line}_1.fastp.fastq.gz -O ${line}_2.fastp.fastq.gz \
    # --detect_adapter_for_pe \
    # -h html/${line}.report.html -j json/${line}.report.json \
    # -q 20 --trim_tail1 1 --trim_tail2 1 -l 20 --thread 4
    echo number of fastp files: `ls /*.fastp.fastq.gz | wc -l`
    echo "\n"
done

echo finished

