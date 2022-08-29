# HPC project - Apriori algorithm implementation using C++ with multi OpenMP multi-threading and SSE.

Environment variables:
export OMP_PROC_BIND="TRUE"

Compile:
g++ apriori.cpp -march=native -O3 -fopenmp
