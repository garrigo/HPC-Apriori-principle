# HPC

Environment variables:
export OMP_PROC_BIND="TRUE"

Compile:
g++ apriori.cpp -march=native -O3 -fopenmp