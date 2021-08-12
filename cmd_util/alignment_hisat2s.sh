#!/bin/sh
# alignment with HISAT2

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check

# below index is from ensemble:
INDEXFILES="/Users/Emma/Documents/Bioinformatics/RefGenome/ENSEMBLE/index/hg38"

# fastq directory
FASTQDIR="/Volumes/HDD24TB/${TARGET}"
# fastp directory
FASTPDIR="/Volumes/HDD24TB/${TARGET}/fastp"
SICKLEDIR="/Volumes/HDD24TB/${TARGET}/trimmed"


# file ids
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021"
FILEIDS="${PROJECTDIR}/${TYPE}/${TARGET}/SRR_Acc_List.txt"


cd ${FASTPDIR}
mkdir -p sam bam bai summary
cat ${FILEIDS} | while read line; do
    echo HISAT2: ${line}
    # hisat2 mapping
    # cf. strand options: https://www.biostars.org/p/258415/
    hisat2 -p 4 --rna-strandness R -x ${INDEXFILES} -U ${line}.fastp.fastq.gz \
    -S sam/${line}.fastp.fastq.sam \
    --summary-file summary/${line}.fastp.hisat2.sum.txt \
    --time
    # hisat2 -p 4 --rna-strandness R -x ${INDEXFILES} -U ${line}.trimmed.fastq.gz \
    # -S sam/${line}.sickle.fastq.sam \
    # --summary-file summary/${line}.sickle.hisat2.sum.txt \
    # --time
    
    # sam to sorted bam
    echo sam to sorted bam...
    samtools sort -@ 4 -O bam -o bam/${line}.fastp.fastq.sort.bam sam/${line}.fastp.fastq.sam
    # samtools sort -@ 4 -O bam -o bam/${line}.sickle.fastq.sort.bam sam/${line}.sickle.fastq.sam
    
    # create index
    echo creating indexed bam file...
    samtools index -@ 4 bam/${line}.fastp.fastq.sort.bam bai/${line}.fastp.fastq.sort.bam.bai
    # samtools index bam/${line}.sickle.fastq.sort.bam 
    echo "\n" 
done

echo finished
