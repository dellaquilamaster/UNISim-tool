#!/bin/bash

MAXEVENTS=100000 #this is the maximum number of events per individual process

if [ -z "$1" ] ; then
  echo "Missing input, nothing done!"
  echo "Please enter number of events"
  exit
fi

TOTEVENTS=$1

EVTPROCESSED=0
ITERATION=0
while [ $EVTPROCESSED -lt $TOTEVENTS ]
do

PROC_NAME=$(printf "UNIS%04d" $ITERATION);
if [ $(( $EVTPROCESSED + $MAXEVENTS )) -le $TOTEVENTS ] ; then
  EVTTOPROCESS=$MAXEVENTS
else
  EVTTOPROCESS=$(( $TOTEVENTS - $EVTPROCESSED ))
fi

cat > submit.sh << EOF
#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -N $PROC_NAME
#$ -e batch_errors/$PROC_NAME.e
#$ -o batch_output/$PROC_NAME.o
#$ -pe mpi 1

module use /share/storage/software/easybuild/modules/all
module load ROOT/6.10.02-foss-2016b-Python-2.7.12

eval './exec_UNISim-tool.exe -events $EVTTOPROCESS -o ${PROC_NAME}_${TOTEVENTS}events.root'

EOF

    qsub submit.sh
#    cat submit.sh
    rm -f submit.sh
    
    sleep 0.2
    
EVTPROCESSED=$(( $EVTPROCESSED + $EVTTOPROCESS))
ITERATION=$(( $ITERATION + 1))

done
