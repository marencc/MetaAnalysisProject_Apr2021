#!bin/sh
# coping command directory to the working directory

# targets
TYPE="inactive"  # check
TARGET="GSE162730"  # check

# directories
PROJECTDIR=/Volumes/HDD24TB/MetaAnalysisProject_Apr2021
COMMANDDIR=${PROJECTDIR}/"command"
WORKINGDIR=${PROJECTDIR}/${TYPE}/${TARGET}/"command"

cp -r ${COMMANDDIR} ${WORKINGDIR}
