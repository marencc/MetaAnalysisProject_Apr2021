#!bin/sh
# scripts for downloading reference genome data

REFDIR="/Users/Emma/Documents/Bioinformatics/RefGenome/"
mkdir -p ${REFDIR}"ENSEMBLE"
ENSEMBLE =${REFDIR}"ENSEMBLE"

cd ${ENSEMBLE}
# release-104 is the latest (Mar 30, 2021 released) as of Aug 2, 2021
lftp ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna

# downloading base assembly fasta file
get Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# gtf file for annotation
cd /pub/release-104/gtf/homo_sapiens
ls
get Homo_sapiens.GRCh38.104.gtf.gz
exit


# decompressing .gz files for future use
pigz -p 4 -d `ls *.gz`
# gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh38.104.gtf.gz

