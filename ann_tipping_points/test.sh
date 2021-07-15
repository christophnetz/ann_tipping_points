#!/bin/bash


P=(1 $(echo 10^0.5 | bc))                             
#R=($(seq 0 0.1 1))
for i in "${P[@]}"
do
#  for j in "${R[@]}"
#do
  echo $i
#  echo $j
#  done
done 
	