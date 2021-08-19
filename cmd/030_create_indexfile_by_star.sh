#!bin/sh
# creating reference genome index files by star

# time
SECONDS=0

# directories
REFDIR=/Volumes/HDD24TB/RefGenome/ENSEMBLE
mkdir -p ${REFDIR}/index/star100bp
STARINDEXDIR=${REFDIR}/index/star100bp

# running STAR
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${STARINDEXDIR} \
--genomeFastaFiles ${REFDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ${REFDIR}/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 100  # default: 100, adjust the number mafter QC

# processed time for one file
h=$(($SECONDS/3600))
m=$((($SECONDS/60)%60))
s=$(($SECONDS%60))
echo processed time: ${h}:${m}:${s}
