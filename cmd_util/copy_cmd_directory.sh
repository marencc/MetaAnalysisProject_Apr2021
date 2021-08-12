#!bin/sh
# coping command directory to the working directory

# targets
TYPE="aging"  # check
TARGET="GSE162730"  # check

# directories
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
mkdir -p ${PROJECTDIR}/${TYPE}/${TARGET}

COMMANDDIR=${PROJECTDIR}/"cmd"
WORKINGDIR=${PROJECTDIR}/${TYPE}/${TARGET}/"cmd"

cp -r ${COMMANDDIR} ${WORKINGDIR}
