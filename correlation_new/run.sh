#!/bin/sh

if [ $# -ne 2 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

#RUN=$1

for RUN in `seq $1 $2`
do

#echo "GEBSort started sorting run $RUN at `date`"
./correlation $RUN
#echo "GEBSort DONE at `date`"
done


#exit
