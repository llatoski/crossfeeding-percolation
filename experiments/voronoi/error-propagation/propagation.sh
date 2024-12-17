#!/bin/bash
# Usage: ./scriptsD TIME DOMAINTYPE
# TIME is the actual time 
# DOMAINTYPE can be ROUGH (non zealots) or SMOOTH (zealots)
# LCFL 29/09/2021
rm lcomponent_frame"$1".dat    
if [ "$1" == 0 ]; then time=9.44  
fi 
if [ "$1" == 1 ]; then time=18.89 
fi 
if [ "$1" == 2 ]; then time=28.33 
fi 
if [ "$1" == 3 ]; then time=37.78 
fi 
if [ "$1" == 4 ]; then time=47.22 
fi 
if [ "$1" == 5 ]; then time=56.67 
fi 
if [ "$1" == 6 ]; then time=66.11 
fi 
if [ "$1" == 7 ]; then time=75.56 
fi 
if [ "$1" == 8 ]; then time=85.00 
fi 
if [ "$1" == 9 ]; then time=94.45 
fi 
if [ "$1" == 10 ]; then time=103.89 
fi 
if [ "$1" == 11 ]; then time=113.34 
fi 
if [ "$1" == 12 ]; then time=122.78 
fi 
if [ "$1" == 13 ]; then time=132.23 
fi 
if [ "$1" == 14 ]; then time=141.67 
fi 
if [ "$1" == 15 ]; then time=151.12 
fi 
if [ "$1" == 16 ]; then time=160.56 
fi 
if [ "$1" == 17 ]; then time=170.01 
fi 
if [ "$1" == 18 ]; then time=179.44 
fi 
if [ "$1" == 19 ]; then time=188.89 
fi 
if [ "$1" == 20 ]; then time=198.33 
fi 
if [ "$1" == 21 ]; then time=207.78 
fi 
if [ "$1" == 22 ]; then time=217.22 
fi 
if [ "$1" == 23 ]; then time=226.67 
fi 
if [ "$1" == 24 ]; then time=236.11 
fi 
if [ "$1" == 25 ]; then time=245.56 
fi 
if [ "$1" == 26 ]; then time=255.00 
fi 
if [ "$1" == 27 ]; then time=264.45 
fi 
if [ "$1" == 28 ]; then time=273.89 
fi 
if [ "$1" == 29 ]; then time=283.34 
fi 
if [ "$1" == 30 ]; then time=292.78 
fi 
if [ "$1" == 31 ]; then time=302.23 
fi 
if [ "$1" == 32 ]; then time=311.67 
fi 
if [ "$1" == 33 ]; then time=321.12 
fi 
if [ "$1" == 34 ]; then time=330.56 
fi 
if [ "$1" == 35 ]; then time=340.01 
fi
count=0; 
for k in $(seq 1 1 1);
do
        output=$(printf 'network.dat')
        cat ../networks/frame"$1"_* > aux
        Nemit=$(cat aux | wc -l)
        Nabsorb=$(awk '{print NF}' aux | sort -nu | tail -n 1)
        awk '$1=='"$time"' {print $2, $3, $4/7, $5/7}' cellstreated.dat >> aux
        gcc errorpropagation_experimental_voronoi_v1.0.c -lm -g -DNEMIT="$Nemit" -DNABSORB="$Nabsorb" -DTHRESHOLD=0.025 -DDT=1 -DSEED="$k"
        ./a.out < aux > $output
        python3 networkstatistics.py >> lcomponent_frame"$1".dat    
        rm -rf network.dat
done
rm aux a.out
