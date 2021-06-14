#include<iostream>
#include<chrono>
#include "apriori3.h"

int main (){
    // int a = (x & 01010101) + (x>>1 & 01010101);
    // int b = (a & 00110011) + (a>>2 & 00110011);
    // int c = (b & 00001111) + (b>>4 & 00001111);
    // std::cout << "OUTPUT: " << c << "\n";
    Apriori apriori;
    apriori.run("chess.dat", 510.0f);
    // apriori.print_single_items();
    // apriori.print_items();

    return 0;
}