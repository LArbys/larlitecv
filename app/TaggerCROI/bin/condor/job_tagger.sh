#!/bin/bash

source /usr/nevis/adm/nevis-init.sh
module load root/06.04.00

let "processid=$1+$2"
sw_dir=$3
cfg_dir=$4
input_dir=$5
out_dir=$6
procspath=$7

scp $procspath procs.txt

let NUM_PROCS=`cat procs.txt | wc -l`
echo "number of processes: $NUM_PROCS"
if [ "$NUM_PROCS" -lt "$processid" ]; then
    echo "No Procces ID to run."
    return
fi

let "proc_line=${processid}+1"
let jobid=`sed -n ${proc_line}p procs.txt`

sleep $jobid
scp -r $sw_dir/larlitecv/app/TaggerCROI/bin/condor ./tmp
#wait
cd tmp
source setup_env.sh

# copy data
cat ${input_dir}/input_larcv_${jobid}.txt | xargs -i scp {} .
wait
cat ${input_dir}/input_larlite_${jobid}.txt | xargs -i scp {} .
wait

# change the filelists to be local
cat ${input_dir}/input_larcv_${jobid}.txt | awk -F/ '{print $NF}' > input_larcv.txt
cat ${input_dir}/input_larlite_${jobid}.txt | awk -F/ '{print $NF}' > input_larlite.txt

# copy cfg file
scp ${cfg_dir}/tagger_${jobid}.cfg .
wait

echo "run_tagger tagger_${jobid}.cfg"

/usr/bin/time -v ${sw_dir}/larlitecv/app/TaggerCROI/bin/./run_tagger tagger_${jobid}.cfg > /dev/null

scp output_*.root ${out_dir}/

ls -lh *
cat input_larcv.txt
cat input_larlite.txt

cd ..

rm -r ./tmp

echo "COMPLETE"
