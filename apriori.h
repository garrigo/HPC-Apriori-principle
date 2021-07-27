#include <vector>
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


class Apriori {
    public:
    std::vector<std::vector<size_t>> transactions;
    std::vector<std::set<size_t>> itemsets;
    std::vector<std::string> single_items;
    std::unordered_map<std::string, size_t> hashmap;
    std::vector<size_t> occurrencies;
    



    inline void print_single_items (){
        for (size_t i=0; i<single_items.size(); i++)
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

    inline void print_items (std::set<std::string>& set , std::string& buff) {
            for (auto str : set){
                std::cout << str <<"    ";
            }
            std::cout << " - Cardinality: " << static_cast<float>(hashmap[buff])/ static_cast<float>(transactions.size()) <<"\n";
        
    }


    void read_data (const std::string input_file) {
        std::ifstream ifs;
        ifs.open(input_file);
            if(!ifs) {
                std::cout << "Input file could not be opened\n";
                exit(0);
            }
            std::string doc_buffer;
            std::vector<std::vector<std::string>> result;
            size_t current_size=0;
            while(!getline(ifs, doc_buffer).eof()){
                std::vector<size_t> line_buffer;
                std::istringstream iss(doc_buffer);
                std::string token;

                while (std::getline(iss, token, ' ')){
                    if (!isspace(token[0])){
                        bool clear=1;
                        // hashmap[token]++;
                        for (size_t i=0; i<single_items.size(); i++){
                            if (single_items[i]==token){
                                clear=0;
                                line_buffer.push_back(i);
                                if (occurrencies.size() < i)
                                    occurrencies.resize(i+1);
                                occurrencies[i]++;
                                i=single_items.size();
                            }
                        }
                        if (clear) { 
                            single_items.push_back(std::move(token));
                            occurrencies.push_back(1);
                            line_buffer.push_back(single_items.size());
                        }
                    }           
                }
                transactions.push_back(std::move(line_buffer));
            }
       
    }

    void singles_prune(float support){
    }

    void singles_merge(float support){
        
        for (size_t i = 0; i<single_items.size()-1; i++){
            if ((static_cast<float>(occurrencies[i]) / (static_cast<float>(transactions.size()))) >=support)
                for (size_t j = i+1; j<single_items.size(); j++){
                    if ((static_cast<float>(occurrencies[j]) / (static_cast<float>(transactions.size()))) >=support){
                        itemsets.push_back({i,j});
                    }
                }
        }
    }


    void map (size_t k){
        std::cout << "ENTER MAP\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        occurrencies.resize(itemsets.size());
        // hashmap.clear();
        //for every itemset
        for (size_t set = 0; set<itemsets.size(); set++){
            occurrencies[set]=0;
            //for every transaction
            for (auto& tx : transactions){
                unsigned int found = 0, cont=1;
                auto item = std::begin (itemsets[set]);
                //for every item in itemset
                while (cont && item != std::end(itemsets[set]) ){
                    cont = 0;
                    //for every item in transaction
                    for (auto item_tx : tx){
                        if ((*item)==(item_tx)){
                            found++;
                            cont = 1;
                            break;
                        }
                    }
                    item++;
                }
                if (found == k){
                    occurrencies[set]++;
                }
            }
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";

        std::cout << "EXIT MAP\n";
    }

    void prune(float support) {
        // //PRUNE
        // // std::cout << "ENTER PRUNE\n";
        
        // std::cout << "ITEMSETS BEFORE PRUNING: " << itemsets.size() << "\n";
        // std::chrono::time_point<std::chrono::system_clock> start, end;
        // start = std::chrono::system_clock::now();

        // auto itemset = std::begin(itemsets);
        // while (itemset != std::end(itemsets)) {
        //     std::string buff;
        //     auto item = std::begin(*itemset);
        //     while (item != std::end(*itemset)){
        //         buff += single_items[*item];
        //         item++;
        //     }
        //     // print_items(*itemset, buff);
        //     if (static_cast<float>(hashmap[buff])/ static_cast<float>(transactions.size()) < support)
        //         itemset = itemsets.erase(itemset);
        //     else
        //         ++itemset;
        // }
        
        // end = std::chrono::system_clock::now();
        // std::chrono::duration<double> elapsed_seconds = end - start;
        // std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
        
        // std::cout << "ITEMSETS AFTER PRUNING: " << itemsets.size() << "\n";
        // // std::cout << "EXIT PRUNE\n";
    }


    void merge (int k, float support){
        // std::cout << "ENTER MERGE\n";

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        if (!itemsets.empty()){
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<std::set<size_t>> temp;
            //for every itemset, try to unite it with another in the itemsets vector
            for (int itemset_x=0; itemset_x<itemsets.size()-1; itemset_x++){
                if ((static_cast<float>(occurrencies[itemset_x]) / (static_cast<float>(transactions.size()))) >=support){
                    for (int itemset_y=itemset_x+1; itemset_y<itemsets.size(); itemset_y++){
                        if ((static_cast<float>(occurrencies[itemset_y]) / (static_cast<float>(transactions.size()))) >=support){
                            std::set<size_t> merged(itemsets[itemset_x]);
                            merged.insert(itemsets[itemset_y].begin(), itemsets[itemset_y].end());
                            if (merged.size() == k)
                                temp.insert(merged);
                        }

                    }
                }
            }
            itemsets.assign( temp.begin(), temp.end() );
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
    }


    void run (const std::string input_file, float support){
        int k=3;
        read_data(input_file);
        singles_merge(support);
        // print_single_items();
        // print_items();
        while (!itemsets.empty()){     
            std::cout << "ENTER PASS N° " << k << "\n\n";   
            map(k);
            // prune(support);
            merge(k, support);
            std::cout << "EXIT PASS N° " << k << "\n\n";
            k++;
            
        }
    }
};