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


Rmax=2.5
Rstep=0.1
Pmax=0.5
Pstep=0.05                            
for i in $(seq 0 $Rstep $Rmax)
do
  for j in $(seq 0 $Pstep $Pmax)
do
  echo $i
  echo $j
  sbatch run_tip_point_loop.sh $i $j
  done
done 
