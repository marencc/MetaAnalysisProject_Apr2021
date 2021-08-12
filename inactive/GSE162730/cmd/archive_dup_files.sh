#!bin/sh
# archive duplicated files

# target
TYPE="inactive"  # check
TARGET="GSE162730"  # check
TARGETDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021/${TYPE}/${TARGET}

# target directory
RNADIR=${TARGETDIR}/fastq_Rnaseq
RIBODIR=${TARGETDIR}/fastq_Riboseq

# creating archive directory
mkdir -p ${RNADIR}/archive ${RIBODIR}/archive
RNA_ARCHIVE=${RNADIR}/archive
RIBO_ARCHIVE=${RIBODIR}/archive

mkdir -p ${RNA_ARCHIVE}/fastp/html ${RNA_ARCHIVE}/fastp/json ${RNA_ARCHIVE}/fastqc \
${RNA_ARCHIVE}/star/finallog/ ${RNA_ARCHIVE}/star/log ${RNA_ARCHIVE}/star/sjout ${RNA_ARCHIVE}/star/unmapped \
${RNA_ARCHIVE}/star_rRnarm/cbam ${RNA_ARCHIVE}/star_rRnarm/finallog ${RNA_ARCHIVE}/star_rRnarm/log ${RNA_ARCHIVE}/star_rRnarm/sjout 

mkdir -p ${RIBO_ARCHIVE}/fastp/html ${RIBO_ARCHIVE}/fastp/json ${RIBO_ARCHIVE}/fastqc \
${RIBO_ARCHIVE}/star_rRnarm/cbam ${RIBO_ARCHIVE}/star_rRnarm/finallog ${RIBO_ARCHIVE}/star_rRnarm/log ${RIBO_ARCHIVE}/star_rRnarm/sjout 

# create duplicated-id file
# check if the files are sorted in the same order before running
# diff --changed-group-format='%<' --unchanged-group-format='' SRR_Acc_List_rna.txt SRR_Acc_List_rna_unique.txt > SRR_dup_rna.txt
# diff --changed-group-format='%<' --unchanged-group-format='' SRR_Acc_List_ribo.txt SRR_Acc_List_ribo_unique.txt > SRR_dup_ribo.txt

RNA_DUP=${TARGETDIR}/SRR_dup_rna.txt
RIBO_DUP=${TARGETDIR}/SRR_dup_ribo.txt

echo moving Rna-seq files...
cat ${RNA_DUP} | while read line; do
    mv ${RNADIR}/${line}* ${RNA_ARCHIVE}
    mv ${RNADIR}/fastp/${line}* ${RNA_ARCHIVE}/fastp
    mv ${RNADIR}/fastp/html/${line}* ${RNA_ARCHIVE}/fastp/html
    mv ${RNADIR}/fastp/json/${line}* ${RNA_ARCHIVE}/fastp/json
    mv ${RNADIR}/fastqc/${line}* ${RNA_ARCHIVE}/fastqc
    mv ${RNADIR}/star/${line}* ${RNA_ARCHIVE}/star
    mv ${RNADIR}/star/finallog/${line}* ${RNA_ARCHIVE}/star/finallog
    mv ${RNADIR}/star/log/${line}* ${RNA_ARCHIVE}/star/log
    mv ${RNADIR}/star/sjout/${line}* ${RNA_ARCHIVE}/star/sjout
    mv ${RNADIR}/star/unmapped/${line}* ${RNA_ARCHIVE}/star/unmapped
    mv ${RNADIR}/star_rRnarm/${line}* ${RNA_ARCHIVE}/star_rRnarm
    mv ${RNADIR}/star_rRnarm/cbam/${line}* ${RNA_ARCHIVE}/star_rRnarm/cbam
    mv ${RNADIR}/star_rRnarm/finallog/${line}* ${RNA_ARCHIVE}/star_rRnarm/finallog
    mv ${RNADIR}/star_rRnarm/log/${line}* ${RNA_ARCHIVE}/star_rRnarm/log
    mv ${RNADIR}/star_rRnarm/sjout/${line}* ${RNA_ARCHIVE}/star_rRnarm/sjout
done  


echo moving Ribo-seq files...
cat ${RIBO_DUP} | while read line; do
    mv ${RIBODIR}/${line}* ${RIBO_ARCHIVE}
    mv ${RIBODIR}/fastp/${line}* ${RIBO_ARCHIVE}/fastp
    mv ${RIBODIR}/fastp/html/${line}* ${RIBO_ARCHIVE}/fastp/html
    mv ${RIBODIR}/fastp/json/${line}* ${RIBO_ARCHIVE}/fastp/json
    mv ${RIBODIR}/fastqc/${line}* ${RIBO_ARCHIVE}/fastqc
    # mv ${RIBODIR}/star/${line}* ${RIBO_ARCHIVE}/star
    # mv ${RIBODIR}/star/finallog/${line}* ${RIBO_ARCHIVE}/star/finallog
    # mv ${RIBODIR}/star/log/${line}* ${RIBO_ARCHIVE}/star/log
    # mv ${RIBODIR}/star/sjout/${line}* ${RIBO_ARCHIVE}/star/sjout
    # mv ${RIBODIR}/star/unmapped/${line}* ${RIBO_ARCHIVE}/star/unmapped
    mv ${RIBODIR}/star_rRnarm/${line}* ${RIBO_ARCHIVE}/star_rRnarm
    mv ${RIBODIR}/star_rRnarm/cbam/${line}* ${RIBO_ARCHIVE}/star_rRnarm/cbam
    mv ${RIBODIR}/star_rRnarm/finallog/${line}* ${RIBO_ARCHIVE}/star_rRnarm/finallog
    mv ${RIBODIR}/star_rRnarm/log/${line}* ${RIBO_ARCHIVE}/star_rRnarm/log
    mv ${RIBODIR}/star_rRnarm/sjout/${line}* ${RIBO_ARCHIVE}/star_rRnarm/sjout
done  
  
echo finished