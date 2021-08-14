# HPC

Environment variables:
export OMP_DISPLAY_ENV="TRUE"
export OMP_PROC_BIND="TRUE"

Compile:
g++ main.cpp -march=native -O3 -fopenmp -o apriori

Execute with perf:
sudo perf stat -d ./apriori