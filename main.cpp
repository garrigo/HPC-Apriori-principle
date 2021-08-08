#include<iostream>
#include<chrono>
// #include "apriori3.h"
#include "apriori.h"
#include "apriori_parallel.h"
constexpr double support = 0.70;
const char file[]="chess.dat";

int main (){
    Apriori apriori;
    ParallelApriori p_apriori;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // apriori.run(file, support);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_seq = end - start;
    start = std::chrono::system_clock::now();

    p_apriori.run(file, support);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_par = end - start;
    
    std::cout <<  "Sequential time: " << elapsed_seconds_seq.count() << "s\nParallel time: " << elapsed_seconds_par.count() << "s\n";


    return 0;
}