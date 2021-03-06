#!/bin/bash


rseed_suffix=123456
mem="1G"
filename="failrate-capture-trans"

mkdir -p logs
mkdir -p data


for N in 10 20 50 100 200
do
	mkdir -p data/${filename}_$N
	mkdir -p logs/${filename}_$N

	for i in {00..99}
	do
		rseed=${N}${i}${rseed_suffix}
		time="$(python primefactors.py $N):00:00"
		jobname="failrate-capture-trans_${N}_${rseed}"
		log="logs/${filename}_$N/log_$rseed.out"
		sbatch --mem=$mem --time=$time --job-name=$jobname --output=$log graham-job-submit.sh $N $rseed $filename
	done
done
		

