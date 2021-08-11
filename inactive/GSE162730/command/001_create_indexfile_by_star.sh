#!bin/sh
# creating reference genome index files by star

# directories
REFDIR="/Volumes/HDD24TB/RefGenome/ENSEMBLE"
mkdir -p "${REFDIR}/index/star50bp"
STARINDEXDIR="${REFDIR}/index/star50bp"

echo creating indexed reference genome file...
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${STARINDEXDIR} \
--genomeFastaFiles ${REFDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ${REFDIR}/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 49  # adjust the number mafter QC
echo finished
