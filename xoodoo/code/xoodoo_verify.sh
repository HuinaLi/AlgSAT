#!/bin/bash
set -e

ShowUsage() {
    echo "verify all xoodoo trails"
    echo ""
    echo "Usage: "
    echo "bash $0 rounds alg_BASE xoodoo_dc_file threads start_i"
    echo ""
    echo "NOTE:"
    echo "by default, bosphorus is at /home/user/anonymous/bosphorus"
    echo "read_dc.py and extension_trail_mess.py should be in alg_BASE"
    echo ""
    echo "Example:"
    echo "bash $0 3 /home/user/anonymous/bosphorus/xoodoo/3rxoodoo/3rxoomess /home/user/anonymous/bosphorus/xoodoo/3rxoodoo/3rxoomess/dc.txt 1 0"
}

if [ "$4" = "" ] || [ "$1" = "--help" ] || [ "$1" = "help" ]; then
    ShowUsage $0
    exit 1
fi

CUR_TIME=`date +%N`
anf_file=xoodoo_verify${CUR_TIME}
rounds=$1
alg_BASE=$2
dc_file=$3
threads=$4
start_i=$5
BASE=`cd $(dirname "$0");pwd`

# check parameters
if [ "${start_i}" = "" ]; then
    start_i=0
fi

file_exist() {
    if [ ! -e "$1" ]; then
        echo "$1 not exits !!!"
        exit 1
    fi
}
file_exist ${alg_BASE}/extension_trail_mess.py
file_exist ${dc_file}
file_exist ${BASE}/solve.sh

echo "bash $0 $@"
echo "start from solution ${start_i}"

generate_anf() {
    python -u ${alg_BASE}/extension_trail_mess.py -i $1 -f ${dc_file} -r ${rounds} >${alg_BASE}/${anf_file}.anf
}

solve() {
    bash ${BASE}/solve.sh ${alg_BASE} ${anf_file} ${threads}
    # remove redundant files
    rm ${alg_BASE}/${anf_file}*
}

main() {
    i=${start_i}
    while true
    do
        echo "**************************************************************************************************"
        echo "start generate anf of solution $i"
        CUR_TIME=`date +%Y/%m/%d/%H:%M:%S`
        echo "current time is: ${CUR_TIME}"
        generate_anf $i
        echo "**************************************************************************************************"
        echo "start solve anf of solution $i"
        solve
        i=$(expr ${i} + 1)
        echo ""
        echo ""
        echo ""
    done
}

main $@