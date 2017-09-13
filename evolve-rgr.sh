#!/bin/bash


rseed_suffix=123456
mem="1G"
filename="rgr-map"


mkdir -p logs
mkdir -p data

transportdist=40 

for N in 10 50 100
do
	mkdir -p data/${filename}_$N
	mkdir -p logs/${filename}_$N

	for i in {00..99}
	do
		rseed=${N}${i}${rseed_suffix}
		time="02:00:00"
		jobname="rgr-map_${N}_${rseed}"
		log="logs/${filename}_$N/log_$rseed.out"
		sbatch --mem=$mem --time=$time --job-name=$jobname --output=$log graham-job-submit.sh $N $rseed $filename $transportdist
	done
done
		

