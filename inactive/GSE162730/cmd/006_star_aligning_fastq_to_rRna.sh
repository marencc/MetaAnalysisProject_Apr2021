#!bin/sh
# removing contaminating rRNA sequence by STAR
# (aligning ribo-seq fastq files to indexed contaminating reference genome sequence)

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check

# directories
PROJECTDIR="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021"
TARGETDIR="${PROJECTDIR}/${TYPE}/${TARGET}/fastq_Rnaseq"
FASTPDIR="${TARGETDIR}/fastp"
# SICKLEDIR="${TARGETDIR}/fastq/sickle"
mkdir -p ${TARGETDIR}/star_rRnarm ${TARGETDIR}/star_rRnarm/cbam
OUTDIR="${TARGETDIR}/star_rRnarm"

# file ids
FILEIDS="${PROJECTDIR}/${TYPE}/${TARGET}/SRR_Acc_List_rna.txt"  # SRR_Acc_List.txt  qclist.txt
FILES=`cat ${FILEIDS} | wc -l`

# contaminating rRna reference indexed file
CONTAMDIR="/Volumes/HDD24TB/RefGenome/CONTAMrRNA/STARindex"

# time
script_started=`date +%s`

cd ${FASTPDIR}
count=1
### STAR
cat ${FILEIDS} | while read line; do
    SECONDS=0
    echo STAR ${count} /${FILES}: ${line} `date "+%m/%d/%Y %H:%M:%S"`

    echo decompressing fastq.gz...
    STAR \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --runMode alignReads --genomeDir ${CONTAMDIR} \
    --readFilesIn ${line}*.fastq.gz \
    --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 \
    --outFileNamePrefix ${OUTDIR}/${line}.rRna.rm. \
    --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
    --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 \
    --alignEndsType EndToEnd > ${OUTDIR}/cbam/${line}.rrna.comtam.bam
    
    echo compressing ${line}.rRna.rm.fastq...
    mv ${OUTDIR}/${line}.rRna.rm.Unmapped.out.mate1 ${OUTDIR}/${line}.rRna.rm.fastq
    pigz -p 4 ${OUTDIR}/${line}.rRna.rm.fastq
    
    # processed time for one file
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
    
    echo processed files: `ls ${OUTDIR}/*.fastq.gz | wc -l` /${FILES}
    count=$((count+1))
    echo "\n" 
done
                                                       
# moving files
cd ${OUTDIR}
mkdir -p log finallog sjout
echo sorting directory...
mv *.final.out finallog
mv *.Log* log
mv *.SJ.out.tab sjout

echo finished
