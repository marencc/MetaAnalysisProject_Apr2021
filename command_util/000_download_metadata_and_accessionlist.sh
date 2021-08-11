#!bin/sh
# download Metadata and Accession List

PROJECT_ID="PRJNA532382"

# Website info
URL=https://trace.ncbi.nlm.nih.gov


# Accession List
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=${PROJECT_ID}&o=acc_s%3Aa

# SRA Run Info (Metadata table)
wget -O Metadata.csv "${URL}/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=${PROJECT_ID}"

# As a parallel to the above example in the Run Selector,
# wget -O ./SRP001599_info.csv "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= SRP001599"
