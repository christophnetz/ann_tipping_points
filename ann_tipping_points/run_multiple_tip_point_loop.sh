#!/bin/bash
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=mult_tip
#SBATCH --output=multiple_tips.log

module load Qt5
export CC=g++
export CXX=g++
make clean
qmake ann_tipping_points.pro
make 

P=(1 3.16288 10 31.6288 100 316.288 1000 3162.88 10000 31628.8 100000)                             
R= $(seq 0 0.1 1)
for i in "${P[@]}"
do
  for j in "${R[@]}"
do
  echo $i
  echo $j
  sbatch run_tip_point_loop.sh $i $j
  done
done 
