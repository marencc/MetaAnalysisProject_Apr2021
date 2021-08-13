#!/bin/sh
# alignment with STAR
# single end

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check

# directories
PROJECTDIR="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021"
TARGETDIR="${PROJECTDIR}/${TYPE}/${TARGET}"
FASTPDIR="${TARGETDIR}/fastq_Rnaseq/fastp"  # check
FILES=`ls ${FASTPDIR}/*.fastq.gz | wc -l`

# STAR output directory
mkdir -p ${TARGETDIR}/fastq_Rnaseq/star
OUTDIR="${TARGETDIR}/fastq_Rnaseq/star"

# file ids
FILEIDS="${TARGETDIR}/SRR_Acc_List_rna.txt"  # SRR_Acc_List.txt  qclist.txt
FILES=`cat ${FILEIDS} | wc -l`  # number of files
# RIBOSEQALLFILEID="${TARGETDIR}/SRR_Acc_List_ribo.txt"
# RIBOSEQALLFILES=`cat ${ROBPSEQALLFILEID} | wc -l`

# ensemble index
INDEXDIR="/Volumes/HDD24TB/RefGenome/ENSEMBLE/index/star50bp"  # check

# time
script_started=`date +%s`

# STAR
cd ${FASTPDIR}
count=1

cat ${FILEIDS} | while read line; do
    SECONDS=0
    echo STAR ${count} /${FILES}: ${line} `date "+%m/%d/%Y %H:%M:%S"`
    
    echo decompressing fastq.gz...
    ### check last 4 lines
    STAR \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --readFilesIn ${line}*.fastq.gz \
    --genomeDir ${INDEXDIR} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUTDIR}/${line}. \
    --quantMode GeneCounts \
    --outReadsUnmapped Fastx \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 0 \
    --outFilterMismatchNmax 2 
    
    # processed time for one file
    echo finished ${line} `date "+%m/%d/%Y %H:%M:%S"`
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
    
    echo processed files: `ls ${OUTDIR}/*.Aligned.sortedByCoord.out.bam | wc -l` /${FILES}
    count=$((count+1))
    echo "\n"
done

echo STAR finished

# moving files
cd ${OUTDIR} # check directory name
mkdir -p log finallog sjout unmapped
echo sorting directory...
mv *.final.out finallog
mv *.Log* log
mv *.SJ.out.tab sjout
mv *.Unmapped.out.mate* unmapped

echo finished
