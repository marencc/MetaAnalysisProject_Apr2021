#!bin/sh
# creating contaminating rRNA reference indexed file by STAR

CONTAMDIR="/Volumes/HDD24TB/RefGenome/CONTAMrRNA"
STARINDEXDIR="${CONTAMDIR}/STARindex"

echo "creating index file for contaminating rRNA fasta..."
STAR \
--runThreadN 4 --runMode genomeGenerate \
--genomeDir ${STARINDEXDIR} \
--genomeFastaFiles ${CONTAMDIR}/contaminating_rRNA.fa
echo finished
