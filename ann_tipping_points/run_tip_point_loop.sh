#!/bin/bash
# Script to run the simulation
# Peregrine directives:
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=single_tip_sim
#SBATCH --output=tip_%j.log

echo "R: "$1
echo "P: "$2
./ann_tipping_points -R $1 -P $2 
