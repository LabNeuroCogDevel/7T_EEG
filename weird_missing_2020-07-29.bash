#!/usr/bin/env bash
set -eu
trap 'e=$?; [ $e -ne 0 ] && echo "$0 exited in error"' EXIT

#
# these are visits that have remarked data but not raw data!?
#  20200729WF  init
mkmissing -1 '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Prep/remarked/1*_2*' -2 '/Volumes/Hera/Raw/EEG/7TBrainMech/1*_2*' | while read ld8; do
 selld8 list|grep -q $ld8 && echo "# should have $ld8, is in db" && continue
 echo -n "# $ld8 not in db. maybe typo. visit db has:  "
 selld8 list|grep eeg -i |grep ${ld8%%_*} |cut -f1 |tr '\n' , 
 selld8 list|grep eeg -i |grep ${ld8##*_} |cut -f1 |tr '\n' ,
 echo
done
