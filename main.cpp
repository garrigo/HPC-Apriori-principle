#include<iostream>
#include<chrono>
// #include "apriori3.h"
#include "apriori.h"

int main (){
    // int a = (x & 01010101) + (x>>1 & 01010101);
    // int b = (a & 00110011) + (a>>2 & 00110011);
    // int c = (b & 00001111) + (b>>4 & 00001111);
    // std::cout << "OUTPUT: " << c << "\n";
    Apriori apriori;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    apriori.run("chess.dat", 0.8f);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";


    return 0;
}