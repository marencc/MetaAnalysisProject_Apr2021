#!/bin/sh
# downloading fastq files by fasterq-dump using SRR Accession ids
# single-end

# time
script_started=`date +%s`

TYPE="aging"  # check
TARGET="GSE164471"  # check

PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
mkdir -p ${PROJECTDIR}/${TYPE}/${TARGET}
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}

# input information
ALLFILEIDS=${TARGETDIR}/SRR_Acc_List.txt
RESTFILEIDS=${TARGETDIR}/SRR_Acc_List_rest.txt
FILES=`cat ${ALLFILEIDS} | wc -l`
RESTFILES=`cat ${RESTFILEIDS} | wc -l`

# output directory
mkdir -p ${TARGETDIR}/fastq
OUTDIR=${TARGETDIR}/fastq
TEMPDIR=${OUTDIR}/temp

count=1
cat $RESTFILEIDS | while read line; do
    SECONDS=0
    echo `date "+%m/%d/%Y %H:%M:%S"` fasterq-dump ${count} /${RESTFILES}: ${line} "(${TARGET})"
    fasterq-dump ${line} --temp ${TEMPDIR} --outdir ${OUTDIR} --threads 8 --progress
    
    echo `date "+%m/%d/%Y %H:%M:%S"` compressing...
    pigz -p 8 ${OUTDIR}/${line}.fastq
    
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
    
    echo processed files: `ls ${OUTDIR}/*.fastq.gz | wc -l` /${ALLFILES}
    count=$((count+1))
    echo "\n"
    done

echo `date "+%m/%d/%Y %H:%M:%S"` fasterq-dump finished "\n"
