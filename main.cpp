#include <iostream>
#include <chrono>
#include "apriori.h"
#include "apriori_sse.h"

int main(int argc, char* argv[])
{
    if(argc!=3)
    {
        std::cout <<"ERROR! Missing/exceeding arguments. Usage: ./apriori <file name> <support>";
        return 1;
    }
    std::string file = argv[1];
    double support = atof(argv[2]);

    SetAprioriSSE set_sse;
    VectorAprioriSSE vector_sse;
    SetApriori set;
    VectorApriori vector;

    std::chrono::time_point<std::chrono::system_clock> start, end;

    //SSE SET SEQUENTIAL
    start = std::chrono::system_clock::now();
    set_sse.run(file, support, false);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> set_sse_seq = end - start;

    set_sse = SetAprioriSSE();

    //SSE SET PARALLEL
    start = std::chrono::system_clock::now();
    set_sse.run(file, support, true);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> set_sse_par = end - start;

    //SSE VECTOR SEQUENTIAL
    start = std::chrono::system_clock::now();
    vector_sse.run(file, support, false);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> vector_sse_seq = end - start;

    vector_sse = VectorAprioriSSE();

    //SSE VECTOR PARALLEL
    start = std::chrono::system_clock::now();
    vector_sse.run(file, support, true);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> vector_sse_par = end - start;

    //NO-SSE SET SEQUENTIAL
    start = std::chrono::system_clock::now();
    set.run(file, support, false);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> set_nosse_seq = end - start;

    set = SetApriori();

    //NO-SSE SET PARALLEL
    start = std::chrono::system_clock::now();
    set.run(file, support, true);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> set_nosse_par = end - start;

    //NO-SSE VECTOR SEQUENTIAL
    start = std::chrono::system_clock::now();
    vector.run(file, support, false);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> vector_nosse_seq = end - start;

    vector = VectorApriori();

    //NO-SSE VECTOR PARALLEL
    start = std::chrono::system_clock::now();
    vector.run(file, support, true);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> vector_nosse_par = end - start;

    std::cout << "\nVERSION                        Execution time (s)\n\n";
    std::cout << "SSE Set Sequential time:       " << set_sse_seq.count() << "\n";
    std::cout << "SSE Set Parallel time:         " << set_sse_par.count() << "\n";
    std::cout << "SSE Vector Sequential time:    " << vector_sse_seq.count() << "\n";
    std::cout << "SSE Vector Parallel time:      " << vector_sse_par.count() << "\n";
    std::cout << "NO-SSE Set Sequential time:    " << set_nosse_seq.count() << "\n";
    std::cout << "NO-SSE Set Parallel time:      " << set_nosse_par.count() << "\n";
    std::cout << "NO-SSE Vector Sequential time: " << vector_nosse_seq.count() << "\n";
    std::cout << "NO-SSE Vector Parallel time:   " << vector_nosse_par.count() << "\n";

    return 0;
}