#!bin/sh
# creating index files by bwa

REFDIR="/Volumes/HDD24TB/RefGenome/ENSEMBLE"
mkdir -p "${REFDIR}/index/bwa"
BWADIR="${REFDIR}/index/bwa"

cd ${REFDIR}
echo creating index file with BWA...
bwa index -p ${BWADIR}/hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa
echo finished

