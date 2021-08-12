#!bin/sh
# scripts for downloading reference genome data

REFDIR="/Users/Emma/Documents/Bioinformatics/RefGenome/"
ENSEMBLE=${REFDIR}"ENSEMBLE/"
mkdir -p ${ENSEMBLE}"index"
ENSEMBLEINDEX=${ENSEMBLE}"index/"

cd ${ENSEMBLE}
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa ${ENSEMBLEINDEX}hg38

