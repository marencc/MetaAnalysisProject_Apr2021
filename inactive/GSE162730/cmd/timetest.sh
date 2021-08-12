#!bin/sh
# scripts for downloading reference genome data

# elapsed time
# timestamp=`date "+%m/%d/%Y %H:%M:%S"`
# echo start: ${timestamp}
# SECONDS=0
# 
# sleep 5
# 
# echo finished: ${timestamp}
# 
# elapsed=${SECONDS}
# h=$(($SECONDS/3600))
# m=$(($SECONDS/60))
# s=$SECONDS
# echo elapsed time: ${h}:${m}:${s}


script_started=`date +%s`
SEQLIBS=(T2 T3 T4 T6 T8 T9)
for seqlib in ${SEQLIBS[@]}; do
    # timestamp=`date "+%m/%d/%Y %H:%M:%S"`
    echo ${seqlib} `date "+%m/%d/%Y %H:%M:%S"`
    SECONDS=0

    sleep 11
    
    # processed time for one file
    echo finished ${line} `date "+%m/%d/%Y %H:%M:%S"`
    h=$(($SECONDS/3600))
    m=$((($SECONDS/60)%60))
    s=$(($SECONDS%60))
    echo processed time: ${h}:${m}:${s}
    
    # elapsed time so far
    elapsed_time=`date +%s`
    elapsed=$((elapsed_time - script_started))
    eh=$(($elapsed/3600))
    em=$((($elapsed/60)%60))
    es=$(($elapsed%60))
    echo elapsed: ${eh}:${em}:${es}
    echo "\n"
done