#include<iostream>
#include<chrono>

#include "apriori.h"
#include "apriori_parallel.h"
#include "apriori_sse.h"

constexpr double support = 0.75;
const char file[]="chess.dat";

int main (){
    SSE_Apriori sse_sequential, sse_parallel;
    Apriori nosse_sequential;
    ParallelApriori nosse_parallel;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //SSE SEQUENTIAL 
    start = std::chrono::system_clock::now();
    sse_sequential.run(file, support, false);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> sse_seq_time= end - start;

    //SSE PARALLEL
    start = std::chrono::system_clock::now();
    sse_parallel.run(file, support, true);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> sse_par_time = end - start;

    //NO-SSE SEQUENTIAL
    start = std::chrono::system_clock::now();
    nosse_sequential.run(file, support);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> nosse_seq_time= end - start;

    //NO-SSE PARALLEL
    start = std::chrono::system_clock::now();
    nosse_parallel.run(file, support);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> nosse_par_time= end - start;

    
    std::cout <<  "SSE: Sequential time: " << sse_seq_time.count() << " - Parallel time: " << sse_par_time.count() << "s\n";
    std::cout <<  "NO-SSE: Sequential time: " << nosse_seq_time.count() << " - Parallel time: " << nosse_par_time.count() << "s\n";

    return 0;
}