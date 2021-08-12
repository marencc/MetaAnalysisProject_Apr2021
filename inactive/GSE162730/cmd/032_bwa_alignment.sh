#!/bin/sh
# alignment with BWA

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check

# directories
PROJECTDIR="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021"
TARGETDIR="${PROJECTDIR}/${TYPE}/${TARGET}"
FASTPDIR="${TARGETDIR}/fastq/fastp"
SICKLEDIR="${TARGETDIR}/fastq/sickle"

# file ids
FILEIDS="${TARGETDIR}/qclist.txt"  # SRR_Acc_List.txt

# ensemble index
INDEXDIR="/Volumes/HDD24TB/RefGenome/ENSEMBLE/index/bwa/hg38"

mkdir -p ${FASTPDIR}/bwa
BWADIR="${FASTPDIR}/bwa"

cd ${BWADIR}
mkdir -p sam bam bai stats
count=1
cat ${FILEIDS} | while read line; do
    echo BWA_${count}: ${line}
    bwa aln ${INDEXDIR} ${FASTPDIR}/${line}.fastp.fastq.gz > sam/${line}.fastp.fastq.sam  # short read
    # bwa bwasw index_name SRR1163656.fastq > SRR1163656.sam  # long read
    # bwa mem index_name SRR1163656.fastq > SRR1163656.sam    # long read + many gaps --> bwa-mem2?

    # sam to sorted bam
    echo sam to sorted bam...
    samtools sort -@ 4 -O bam -o bam/${line}.fastp.fastq.sort.bam sam/${line}.fastp.fastq.sam
    # samtools sort -@ 4 -O bam -o bam/${line}.sickle.fastq.sort.bam sam/${line}.sickle.fastq.sam
    
    # create index
    echo creating indexed bam file...
    samtools index -@ 4 bam/${line}.fastp.fastq.sort.bam bai/${line}.fastp.fastq.sort.bam.bai
    # samtools index bam/${line}.sickle.fastq.sort.bam 
    
    # alignment stats
    samtools flagstat bam/${line}.fastp.fastq.sort.bam stats/${line}.fastp.fastq.sort.bam.stats
    
    count=$((count+1))
    echo "\n" 
done

echo finished
