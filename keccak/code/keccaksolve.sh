#!/bin/bash

set -e

ShowUsage() {
    echo "Usage: "
    echo "bash $0 alg_BASE anf_filename threads"
    echo ""
    echo "NOTE:"
    echo "by default, bosphorus is at /home/user/anonymous/bosphorus"
    echo "anf_file should be in alg_BASE"
    echo ""
    echo "Example: (gimli)"
    echo "bash $0 /home/user/anonymous/bosphorus/gimli 8rgimli 20"
}

if [ "$3" = "" ] || [ "$1" = "--help" ] || [ "$1" = "help" ]; then
    ShowUsage $0
    exit 1
fi

bos_BASE=/home/user/anonymous/bosphorus
alg_BASE=$1
filename=$2
threads=$3
cd ${bos_BASE}
echo "bash $0 $1 $2 $3"

# create anf

echo "start creating anf"
nohup python -u /home/user/anonymous/bosphorus/keccak/4rkeccak.py -f /home/user/anonymous/bosphorus/keccak/trails_1600.txt </dev/null >${alg_BASE}/${filename}.anf 2>&1
echo "output anf file: ${alg_BASE}/${filename}.anf"
echo "end creating anf"


# start time log
START_TIME=`date +%Y/%m/%d/%H:%M:%S`
ST=`date +%s.%N`
echo "solve start time is: ${START_TIME}"

# solve
${bos_BASE}/build/bosphorus --anfread ${alg_BASE}/${filename}.anf --anfwrite ${alg_BASE}/${filename}_out.anf \
            --cnfwrite ${alg_BASE}/${filename}.cnf --solvewrite ${alg_BASE}/${filename}_solution -t ${threads}

# end time log
ED=`date +%s.%N`
END_TIME=`date +%Y/%m/%d/%H:%M:%S`
echo "solve end time is: ${END_TIME}"
set +e
EXECUTING_TIME=$(printf "%.6f" `echo "${ED} - ${ST}" | bc`)
set -e
echo "solve exec time is: ${EXECUTING_TIME}s"