#/bin/bash
DATE=`date +%F`
old="$IFS"
IFS="-"
args_str="$*"
IFS=$old
FILENAME="cachegrind.out.${DATE}_${args_str}_%n"
nice -n +14 valgrind --tool=callgrind --callgrind-out-file=$FILENAME C/csteen $*