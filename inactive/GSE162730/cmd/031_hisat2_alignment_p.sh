#!/bin/sh
# alignment

# targets
TYPE="inactive/"  # check
TARGET="GSE130722/"  # check

# index file is obtained from http://daehwankimlab.github.io/hisat2/download/#h-sapiens
# UCSC hg38
# https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
INDEXDIR="/Users/Emma/Documents/Bioinformatics/RefGenome/UCSC/RefHg38/hg38/"

# fastqc directory
HDD="/Volumes/HDD24TB/"
FASTQDIR=${HDD}${TARGET}

# file ids
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/"
FILEIDS=${PROJECTDIR}${TYPE}${TARGET}"SRR_Acc_List.txt"

cd $FASTQDIR
cat $FASTQIDS | while read line; do
    echo HISAT2: ${line}
    #hisat2 mapping
    hisat2 -x $INDEXDIR/genome -1 ${line}_1.fastq -2 ${line}_2.fastq -p 8 -S ${line}.sam \
    --time --summary-file ${line}_sum.txt
    #sam to sorted bam
    echo sam to sorted bam...
    samtools sort -@ 4 -O bam -o ${line}.sort.bam ${line}.sam
    #create index
    echo creating indexed bam file...
    samtools index ${line}.sort.bam
    echo "\n"
done

### Moving files
echo moving files...
mkdir -p sam bam bai summary
mv *.sam sam
mv *.bam bam
mv *.bai bai
mv *_sum.txt summary

echo finished
