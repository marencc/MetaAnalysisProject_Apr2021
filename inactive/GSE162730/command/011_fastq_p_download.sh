#!/bin/sh
### Obtaining fastq files by fasterq-dump using SRR Accession ids.

HDD="/Volumes/HDD24TB/"
TARGET="GSE113165"  # check
mkdir -p ${HDD}${TARGET}

PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/"
TYPE="inactive"  # check
FILEIDS=${PROJECTDIR}${TYPE}"/"${TARGET}"/SRR_Acc_List.txt"
OUTDIR=${HDD}${TARGET}
TEMPDIR=${OUTDIR}"/temp"

cat $FILEIDS | while read line; do
    echo Downloading: ${line}
    fasterq-dump ${line} --split-files --split-3 \
    --temp ${TEMPDIR} --outdir ${OUTDIR} --threads 4 --progress
    pigz -p 4 ${OUTDIR}"/"${line}*fastq
    echo "\n";
done
echo finished
