#!bin/sh

# change directory: reference genome
cd /Users/Emma/Documents/Bioinformatics/RefGenome/ENSEMBLE/

# converting gtf file to bed file
gtf2bed < Homo_sapiens.GRCh38.104.gtf > Homo_sapiens.GRCh38.104.bed
