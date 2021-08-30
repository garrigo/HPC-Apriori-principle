#include <iostream>
#include <chrono>
#include "apriori.h"
#include "apriori_sse.h"

int main(int argc, char* argv[])
{
    if(argc!=4)
    {
        std::cout <<"ERROR! Missing/exceeding arguments. Usage: ./apriori <file name> <support>";
        return 1;
    }
    std::string file = argv[1];
    const double support = atof(argv[2]);
    const int max_threads = std::max(1, std::min(atoi(argv[3]), omp_get_max_threads()));

    AprioriSSE sse_apriori;
    SyncAprioriSSE sync_sse;
    Apriori apriori;
    SyncApriori sync_apriori;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    std::cout << "EXECUTING WITH " << max_threads << " THREADS\n"; 
    std::cout << "\nVERSION                        Execution time (s)\n\n";

    //SSE SEQUENTIAL
    start = std::chrono::system_clock::now();
    sse_apriori.run(file, support, 1);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "SSE Sequential time:           " << time.count() << "\n";

    sse_apriori = AprioriSSE();

    //SSE PARALLEL
    start = std::chrono::system_clock::now();
    sse_apriori.run(file, support, max_threads);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "SSE Parallel time:             " << time.count() << "\n";
    

    //SSE SYNC PARALLEL
    start = std::chrono::system_clock::now();
    sync_sse.run(file, support, max_threads);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "SSE Sync Parallel time:        " << time.count() << "\n";

    //NO-SSE SEQUENTIAL
    start = std::chrono::system_clock::now();
    apriori.run(file, support, 1);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "NO-SSE Sequential time:        " << time.count() << "\n";

    apriori = Apriori();

    //NO-SSE PARALLEL
    start = std::chrono::system_clock::now();
    apriori.run(file, support, max_threads);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "NO-SSE Parallel time:          " << time.count() << "\n";



    //NO-SSE SYNC PARALLEL
    start = std::chrono::system_clock::now();
    sync_apriori.run(file, support, max_threads);
    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "NO-SSE Sync Parallel time:     " << time.count() << "\n";

    return 0;
}