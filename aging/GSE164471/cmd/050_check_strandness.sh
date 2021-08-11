#!bin/sh
# checking strandness

# targets
TYPE="inactive"  # check
TARGET="GSE113165"  # check
ID="SRR7007949.sort.bam"  # check

# input: bam file directory 
BAMDIR="/Volumes/HDD24TB/${TARGET}/re/bam/"

# ref: bed file directory
BEDFILE="/Users/Emma/Documents/Bioinformatics/RefGenome/ENSEMBLE/Homo_sapiens.GRCh38.104.bed"

# outdir
OUTDIR="/Users/Emma/Documents/Bioinformatics/DEG/MetaAnalysisProject_Apr2021/${TYPE}/${TARGET}/"

# for bam in `ls ${FASTQDIR}` | while read line; do
infer_experiment.py -r ${BEDFILE} -i ${BAMDIR}${ID} > ${OUTDIR}${ID}_strand_stat.txt

