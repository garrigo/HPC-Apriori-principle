#include <vector>
#include <deque>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <chrono>
#include <bitset>
#include <omp.h>



using IndexType = unsigned int;

struct VectorHash {
    inline IndexType operator()(const std::vector<IndexType>& v) const {
        std::hash<IndexType> hasher;
        IndexType seed = 0;
        for (IndexType i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

using U_VectorSet = std::unordered_set<std::vector<IndexType>, VectorHash>;


class Apriori {
    public:
    std::vector<std::vector<IndexType>> transactions;
    std::vector<std::vector<IndexType>> itemsets;
    std::vector<std::string> single_items;
    std::vector<IndexType> occurrencies;
    

    inline void print_single_items (){
        for (IndexType i=0; i<single_items.size(); i++)
            std::cout <<"Element: " << single_items[i] << " - Cardinality: "<< occurrencies[i] << "\n";
        std::cout << "SINGLE ITEMS Total size: " << single_items.size() << "\n";
    }

    inline void print_items () {
        for (auto set : itemsets){
            for (auto index : set){
                std::cout << single_items[index] <<"    ";
            }
            std::cout << "\n";
        }
        std::cout << "ITEMSETS total size: " << itemsets.size() << "\n";
    }


    void read_data (const std::string input_file, bool parallel=1) {
        std::ifstream ifs;
        ifs.open(input_file);
            if(!ifs) {
                std::cout << "Input file could not be opened\n";
                exit(0);
            }
            std::string doc_buffer;
            std::vector<std::vector<std::string>> result;
            IndexType current_size=0;
            while(!getline(ifs, doc_buffer).eof()){
                std::vector<IndexType> line_buffer;
                std::istringstream iss(doc_buffer);
                std::string token;

                while (std::getline(iss, token, ' ')){
                    if (!isspace(token[0])){
                        bool clear=1;
                        for (IndexType i=0; i<single_items.size(); i++){
                            if (single_items[i]==token){
                                clear=0;
                                line_buffer.push_back(i);
                                if (occurrencies.size() < i)
                                    occurrencies.resize(i+1);
                                occurrencies[i]++;
                                i=single_items.size();
                                break;
                            }
                        }
                        if (clear) { 
                            single_items.push_back(std::move(token));
                            occurrencies.push_back(1);
                            line_buffer.push_back(single_items.size());
                        }
                    }           
                }
                std::sort(line_buffer.begin(), line_buffer.end());
                transactions.push_back(std::move(line_buffer));
            }
       
    }


    void singles_merge(double support, bool parallel=1){
        IndexType size = transactions.size();
        #pragma parallel for
        for (IndexType i = 0; i<single_items.size()-1; i++){
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >=support)
                // #pragma omp for schedule(static) nowait
                for (IndexType j = i+1; j<single_items.size(); j++){
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >=support){
                        itemsets.push_back({i,j});
                    }
                }
            // if (i%1000==0){
            //     #pragma omp barrier
            // }
        }
    }


    void map (IndexType k, bool parallel=1){
        // std::cout << "ENTER MAP\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        occurrencies.resize(itemsets.size());
        //for every itemset
        #pragma omp parallel for num_threads(4)
        for (IndexType set=0; set<itemsets.size(); set++){
        // for (const auto set : itemsets){
            occurrencies[set]=0;
            //for every transaction
            for (IndexType tx=0; tx<transactions.size(); tx++){
                if (transactions[tx].size()>=k){
                    IndexType found = 0, cont=1, tx_cursor=0;
                    auto item = itemsets[set].begin();
                    //for every item in itemset
                    while (cont && item != itemsets[set].end() ){
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor<transactions[tx].size()){
                            // if ((*item)<(tx[tx_cursor])){
                            //     tx_cursor=tx.size();
                            // }
                            // else 
                            if ((*item)==(transactions[tx][tx_cursor])){
                                ++found;
                                cont = 1;
                            }

                            ++tx_cursor;
                        }
                        ++item;
                    }
                    if (found == itemsets[set].size()){
                        occurrencies[set]++;
                    }
                }
            }
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Map time: " << elapsed_seconds.count() << "s\n";

        // std::cout << "EXIT MAP\n";
    }

    void prune (float support){
        // std::cout << "ITEMSETS BEFORE PRUNING: " << itemsets.size() << "\n";
        // std::chrono::time_point<std::chrono::system_clock> start, end;
        // start = std::chrono::system_clock::now();        
        // IndexType occ = 0;
        // IndexType size = transactions.size();
        // auto item_set = std::begin(itemsets);
        // while (item_set != std::end(itemsets)){
        //     if ((static_cast<float>(occurrencies[occ]) / (static_cast<float>(size))) < support){
        //         item_set = itemsets.erase(item_set);
        //     }
        //     else
        //         ++item_set;
        //     ++occ;
        // }
        // end = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds = end - start;
        // std::cout <<  "Prune time: " << elapsed_seconds.count() << "s\n";
        // std::cout << "ITEMSETS AFTER PRUNING: " << itemsets.size() << "\n";
    }

    void merge (unsigned int k, double support){
        // std::cout << "ENTER MERGE\n";

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        if (!itemsets.empty()){
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            U_VectorSet temp;
            IndexType size = transactions.size();
            //for every itemset, try to unite it with another in the itemsets vector
            // #pragma omp parallel num_threads(4)

            // std::cout << "OCCURRENCIES SIZE: " << occurrencies.size() << "\n";
            #pragma omp parallel num_threads(4)
            for (IndexType i=0; i<itemsets.size()-1; i++){
                if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                {
                    #pragma omp for schedule(static) nowait
                        for(IndexType j=i+1; j<itemsets.size(); j++){    
                            if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                            {   
                                std::vector<IndexType> merged(k);
                                IndexType m=0, v1=0, v2=0, distance=0;
                                while (distance<3 && m<k){
                                    if (v2==k-1 || (v1!=k-1 && itemsets[i][v1] < itemsets[j][v2])){
                                        merged[m] = itemsets[i][v1];
                                        ++v1;
                                        ++distance;
                                    }
                                    else if (v1==k-1 || (v2!=k-1 && itemsets[i][v1] > itemsets[j][v2])){
                                        merged[m] = itemsets[j][v2];
                                        ++v2;
                                        ++distance;
                                    }
                                    else {
                                        merged[m] = itemsets[j][v2];
                                        ++v1;
                                        ++v2;
                                    }
                                    ++m;
                                }
                                if (distance==2){
                                    #pragma omp critical (write_temp)
                                    {temp.insert(merged);}
                                }
                            }
                        }
                    if (i%100 == 0){
                        #pragma omp barrier
                    }
                }
            }
            itemsets.resize(temp.size()); int i=0;
            for(auto& element : temp){ itemsets[i] = element; ++i;}
            std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Merge time: " << elapsed_seconds.count() << "s\n";
    }


    void run (const std::string & input_file, double support){
        unsigned int k=2;
        read_data(input_file);
        singles_merge(support);
        // print_single_items();
        // print_items();
        while (!itemsets.empty()){     
            std::cout << "ENTER PASS N° " << k << "\n";   
            map(k);
            prune(support);
            ++k;
            merge(k, support);
            std::cout << "EXIT PASS N° " << k-1 << "\n\n";
        }
    }
};