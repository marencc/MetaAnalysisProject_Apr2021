#!/bin/sh
# trimming low quality reads by fastp
# adapter trimming is enabled by default
# detecting adapter is TRUE for single end by default

# time
script_started=`date +%s`

# targets
TYPE="aging"  # check
TARGET="GSE164471"  # check
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}

# fastq directory
FASTQDIR=${TARGETDIR}/fastq

# output directory
mkdir -p ${TARGETDIR}/fastp
FASTPDIR=${TARGETDIR}/fastp

# file ids
FILEIDS=${TARGETDIR}/SRR_Acc_List.txt
FILES=`cat ${FILEIDS} | wc -l`

cd ${FASTPDIR}
mkdir -p html json
count=1
cat $FILEIDS | while read line; do
    SECONDS=0
    echo `date "+%m/%d/%Y %H:%M:%S"` fastp ${count} /${FILES}: ${line} "(${TARGET})"
    # single end
    # -3: enable per read cutting by quality in tail (3'), default is disabled
    # (WARNING: this will interfere deduplication for SE data)
    fastp -i ${FASTQDIR}/${line}.fastq.gz -o ${line}.fastp.fastq.gz \
    -h html/${line}.report.html -j json/${line}.report.json --thread 6 -q 20 
    
    # pair end
    # fastp -i ${FASTQDIR}/${line}_1.fastq.gz -I ${FASTQDIR}/${line}_2.fastq.gz -3 \
    # -o ${line}_1.fastp.fastq.gz -O ${line}_2.fastp.fastq.gz \
    # --detect_adapter_for_pe \
    # -h html/${line}.report.html -j json/${line}.report.json \
    # -q 20 --trim_tail1 1 --trim_tail2 1 -l 20 --thread 4
    echo `date "+%m/%d/%Y %H:%M:%S"` ${line} finished
    
    # processed time for one file
    echo `date "+%m/%d/%Y %H:%M:%S"` finished ${line}
    h=$(($SECONDS/3600))
    m=$((($SECONDS/60)%60))
    s=$(($SECONDS%60))
    echo processed time: ${h}:${m}:${s}
    
    # elapsed time so far
    elapsed_time=`date +%s`
    elapsed=$((elapsed_time - script_started))
    eh=$(($elapsed/3600))
    em=$((($elapsed/60)%60))
    es=$(($elapsed%60))
    echo elapsed: ${eh}:${em}:${es}
    
    echo number of fastp files: `ls *.fastp.fastq.gz | wc -l` /${FILES}
    count=$((count+1))
    echo "\n"
done

echo `date "+%m/%d/%Y %H:%M:%S"` fastp finished "\n"

