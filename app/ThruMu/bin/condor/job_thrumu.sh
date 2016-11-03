#!/bin/bash

source /usr/nevis/adm/nevis-init.sh
module load root/06.04.00

let "jobid=$1+$2"
sw_dir=$3
cfg_dir=$4
input_dir=$5
out_dir=$6

sleep $jobid
scp -r $sw_dir/larlitecv/app/ThruMu/bin/condor ./tmp
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
scp ${cfg_dir}/bmt_${jobid}.cfg .
wait

echo "thrumu bmt_${jobid}.cfg input_larcv.txt input_larlite.txt"

/usr/bin/time -v thrumu bmt_${jobid}.cfg input_larcv.txt input_larlite.txt
#thrumu bmt_${jobid}.cfg input_larcv.txt input_larlite.txt
#wait

scp output_*.root ${out_dir}/
#wait

ls -lh *
cat input_larcv.txt
cat input_larlite.txt

cd ..

#rm -r ./tmp

echo "COMPLETE"
