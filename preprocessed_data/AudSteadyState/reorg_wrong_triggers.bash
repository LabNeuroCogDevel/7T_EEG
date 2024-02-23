#!/usr/bin/env bash
# DONT USE
#
# some recordings in ICAwhole have unexpected triggers
# 
# these have been corrected by remark files
# and we want those to be the input to future processing
# one solution is to move the originals to a "wrong" file name
# and move remarked as if nothing ever happend
# ... a better way might be to modify the processing code to look for 
# icpru_remarked.set first and only use icpru.set only if it doesn't exist
# then we dont need this script
#
#
# 20221003WF - init
# 20221004WF - DONT USE. matlab checks for remarked file instead of needing to rename 
#
[ -v DRYRUN ] && DRYRUN=echo || DRYRUN=
warn(){ echo "$@" >&2; }

id_all() { ls 1*_2*|grep -oP '\d{5}_\d{8}'|sort -u; }
id_four_files() { ls 1*_2*.set|grep -oP '\d{5}_\d{8}'|sort|uniq -c | awk '($1==4){print $2}' ; }

rename_files(){
   local id="$1"

   # if there's no remarked we dont want to do anything
   test ! -r "${id}_ss_Rem_rerefwhole_ICA_pesos_remarked.set" &&
      echo "# skip: no '$_'. already renamed?" && return 0

   # move originals to bad
   rename -n 's/.(set|fdt)/_wrongTriggers.$1/' $(ls "$id"*|grep -v remarked)

   # TODO rename to _.$1 or something better instead
   # OR keep _remarked
   rename -n 's/_remarked.(set|fdt)/.$1/' "$id"*remarked.*
}

rename_all() {
   cd /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/AudSteadyState/AfterWhole/ICAwholeClean
 for id in $(id_four_files); do  
    rename_files "$id"
 done
 return 0
}

# if not sourced (testing), run as command
if ! [[ "$(caller)" != "0 "* ]]; then
  echo "no move need. matlab looks for remarked" && exit 1
  set -euo pipefail
  trap 'e=$?; [ $e -ne 0 ] && echo "$0 exited in error $e"' EXIT
  rename_all "$@"
  exit $?
fi

####
# testing with bats. use like
#   bats ./reorg_wrong_triggers.bash --verbose-run
####
function fourfiles_test { #@test 

   cd /Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/AudSteadyState/AfterWhole/ICAwholeClean
   run id_four_files
   [[ $output =~ 11875_20220721 ]]

   [ 0 -eq 1 ]
}
