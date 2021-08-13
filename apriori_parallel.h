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
constexpr unsigned int byte_size=sizeof(IndexType);


class ParallelApriori {
    std::vector<std::vector<IndexType>> transactions;
    std::vector<std::vector<IndexType>> itemsets;
    std::vector<std::string> single_items;
    std::vector<IndexType> occurrencies;

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
    using VectorSet = std::unordered_set<std::vector<IndexType>, VectorHash>;

    inline void print_single_items (){
        for (IndexType i=0; i<single_items.size(); i++){
            //std::cout <<"Element: " << single_items[i] << " - Cardinality: "<< occurrencies[i] << "\n";
        }
        //std::cout << "SINGLE ITEMS Total size: " << single_items.size() << "\n";
    }

    inline void print_items () {
        for (auto set : itemsets){
            for (auto index : set){
                //std::cout << single_items[index] <<"    ";
            }
            //std::cout << "\n";
        }
        //std::cout << "ITEMSETS total size: " << itemsets.size() << "\n";
    }


    void read_data (const std::string input_file) {
        std::ifstream ifs;
        ifs.open(input_file);
            if(!ifs) {
                //std::cout << "Input file could not be opened\n";
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


    void singles_merge(double support){
        const IndexType tx_size = transactions.size(), single_size = single_items.size();
        // const IndexType cache_regulator = std::max((500000/byte_size), 1U);
        #pragma parallel for schedule(dynamic,1)
        for (IndexType i = 0; i<single_size-1; i++){
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(tx_size))) >=support)
                // #pragma omp for schedule(static) nowait
                for (IndexType j = i+1; j<single_size; j++){
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(tx_size))) >=support){
                        itemsets.push_back({i,j});
                    }
                }
            // if (i%cache_regulator==0){
            //     #pragma omp barrier
            // }
        }
        std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
    }

    void map1 (IndexType k){
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        const IndexType cache_regulator = std::max((200000/byte_size*(unsigned int)transactions[0].size()), (unsigned int)1);
        occurrencies.resize(itemsets.size());
        // occurrencies = std::vector<IndexType>(itemsets.size(), 0);
        #pragma omp parallel for
        for (IndexType i=0; i<occurrencies.size(); i++)
            occurrencies[i]=0;
        //for every itemset
        #pragma omp parallel 
        
        for (IndexType tx=0; tx<transactions.size(); tx++){
        // for (IndexType set=0; set<itemsets.size(); set++){
        // for (const auto set : itemsets){
            // #pragma omp single
            // {occurrencies[set]=0;}
            //for every transaction
            #pragma omp for schedule(static) nowait
            for (IndexType set=0; set<itemsets.size(); set++){
            // for (IndexType tx=0; tx<transactions.size(); tx++){
                // if (transactions[tx].size()>=k){
                    IndexType found = 0, cont=1, tx_cursor=0;
                    auto item = itemsets[set].begin();
                    //for every item in itemset
                    while (cont && item != itemsets[set].end() ){
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor<transactions[tx].size()){
                            if ((*item)<(transactions[tx][tx_cursor])){
                                break;
                            }
                            else 
                            if ((*item)==(transactions[tx][tx_cursor])){
                                ++found;
                                cont = 1;
                            }

                            ++tx_cursor;
                        }
                        ++item;
                    }
                    if (found == itemsets[set].size()){
                        #pragma omp atomic
                        occurrencies[set]++;
                    }
                // }
            }
            if (tx%cache_regulator==0){
                #pragma omp barrier
            }
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        //std::cout <<  "Map time: " << elapsed_seconds.count() << "s\n";
    }
    
    void map (IndexType k){
        // std::chrono::time_point<std::chrono::system_clock> start, end;
        // start = std::chrono::system_clock::now();

        // const IndexType cache_regulator = std::max((500000/byte_size), (unsigned int)1);
        occurrencies.resize(itemsets.size());
        //for every itemset
        #pragma omp parallel for
        for (IndexType set=0; set<itemsets.size(); set++){
            occurrencies[set]=0;
            //for every transaction
            for (IndexType tx=0; tx<transactions.size(); tx++){
                if (transactions[tx].size()>=k){
                    IndexType found = 0, cont=1, tx_cursor=0;
                    auto item = itemsets[set].begin();
                    //for every item in itemset
                    while (cont && item != itemsets[set].end() ){
                        cont = 0;
                        //for every item in a transaction
                        while (!cont && tx_cursor<transactions[tx].size()){
                            if ((*item)<(transactions[tx][tx_cursor])){
                                break;
                            }
                            else 
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

        // end = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds = end - start;
        //std::cout <<  "Map time: " << elapsed_seconds.count() << "s\n";

    }


    void merge (unsigned int k, double support){
        // std::chrono::time_point<std::chrono::system_clock> start, end;
        // start = std::chrono::system_clock::now();

        if (!itemsets.empty()){
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<std::vector<IndexType>> temp;
            std::vector<std::vector<IndexType>> v_temp;
            IndexType size = transactions.size();
            IndexType itemsets_size = itemsets.size();
            // const IndexType cache_regulator = std::max((500000/byte_size*(k-1)), (unsigned int)1);
            //for every itemset, try to unite it with another in the itemsets vector

            // #pragma omp parallel
            #pragma omp parallel for schedule(dynamic)
            for (IndexType i=0; i<itemsets_size-1; i++){
                if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                {
                    // #pragma omp for schedule(static) nowait
                        for(IndexType j=i+1; j<itemsets_size; j++){    
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
                                    #pragma omp critical (merge_write)
                                    {
                                        if (temp.insert(merged).second)
                                            v_temp.push_back(merged);
                                    }
                                }
                            }
                        }
                    // if (i%cache_regulator == 0){
                    //     #pragma omp barrier
                    // }
                }
            }
            itemsets.swap(v_temp);
            // itemsets.resize(temp.size(), std::vector<IndexType>(k));
            // IndexType i=0;
            // for(auto& element : temp){
            //     // #pragma omp task firstprivate(i)
            //     {itemsets[i] = element;}
            //     ++i;
            // }
            std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }

        // end = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds = end - start;
        //std::cout <<  "Merge time: " << elapsed_seconds.count() << "s\n";
    }

public:
    void run (const std::string & input_file, double support){
        unsigned int k=2;
        read_data(input_file);
        singles_merge(support);
        while (!itemsets.empty()){      
            map(k);
            ++k;
            merge(k, support);
        }
    }
};