#!bin/sh
# download Metadata and Accession List

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check
PROJECT_ID="PRJNA682757"

# output directory
# PROJECTDIR="/Volumes/HDD24TB/MetaAnalysisProject_Apr2021"
PROJECTDIR="/Users/Emma/Documents/Bioinformatics/MetaAnalysisProject_Apr2021"
mkdir -p ${PROJECTDIR}/${TYPE}/${TARGET}
TARGETDIR=${PROJECTDIR}/${TYPE}/${TARGET}

# Website info
URL=https://trace.ncbi.nlm.nih.gov

# SRA Run Info (Metadata table)
wget -O ${TARGETDIR}/metadata.csv "${URL}/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=${PROJECT_ID}"

# As a parallel to the above example in the Run Selector,
# wget -O ./SRP001599_info.csv "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= SRP001599"
