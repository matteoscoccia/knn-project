#!/bin/bash


mpicc -std=c99 -o knn_parallel parallel.c -lm


potenzaMassimaProcessi=6
potenzaMassimaK=16
massimoValoreK=20


for ((i=1; i<=potenzaMassimaProcessi; i++))
do
    p=$((2**i)) 
    for ((j=1; j<=potenzaMassimaK; j++))
    do
        n=$((2**j)) 
        for ((h=5; h<=massimoValoreK; h+=5))
        do
            k=$((h))
            mpirun ./knn_parallel $p $n $k 
        done
    done
done
