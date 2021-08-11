#!bin/sh
# sorting RNA-seq and ribo-seq libraries

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check

# fastq directory
PROJECTDIR="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021"
TARGETDIR="${PROJECTDIR}/${TYPE}/${TARGET}/fastqc"

cd ${TARGETDIR}
echo creating directories...
mkdir -p RNAseq Riboseq

RNASEQIDS="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE162730/SRR_Acc_List_rna.txt"
RIBOSEQIDS="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/inactive/GSE162730/SRR_Acc_List_ribo.txt"

echo moving rna-seq files...
cat ${RNASEQIDS} | while read line; do
    mv ${line}* RNAseq
done

echo moving ribo-seq files...
cat ${RIBOSEQIDS} | while read line; do
    mv ${line}* RIBOseq
done

echo finished
