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

P=(1 $((10**0.5)) 10 $((100**0.5)) 100 $((1000**0.5)) 1000 $((10000**0.5)) 10000 $((100000**0.5)) 100000)                             
R= $(seq 0 0.1 1)
for i in "${P[@]}"
do
  for j in "${R[@]}"
do
  echo $i
  echo $j
#  sbatch run_tip_point_loop.sh $i $j
  done
done 
