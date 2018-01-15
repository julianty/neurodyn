#!/bin/bash
#$ -t 1-150
#$ -N Rea_820
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -m beas
#$ -o ./output
#$ -e ./error
#$ -q batch.q
./simple_nakl_cpp $SGE_TASK_ID
