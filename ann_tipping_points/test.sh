#!/bin/bash


P=(1 10**0.5 10 100**0.5 100 1000**0.5 1000 10000**0.5 10000 100000**0.5 100000)                             
R= $(seq 0 0.1 1)
for i in "${P[@]}"
do
  for j in "${R[@]}"
do
  echo $i
  echo $j
  done
done 
	