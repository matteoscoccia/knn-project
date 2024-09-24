#!/bin/bash


mpicc -std=c99 -o knn_sequential sequential.c -lm


potenzaMassimaN=16
massimoValoreK=20

for ((j=1; j<=potenzaMassimaN; j++))
do
    n=$((2**j)) 
    for ((h=5; h<=massimoValoreK; h+=5))
    do
        k=$((h))
        ./knn_sequential $n $k
    done
done
