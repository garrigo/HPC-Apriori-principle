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
    std::set<std::vector<size_t>> itemsets;
    std::vector<std::string> single_items;
    std::unordered_map<std::string, size_t> hashmap;
    std::vector<size_t> occurrencies;
    

    int i = sizeof(unsigned int);

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
                        itemsets.insert({i,j});
                    }
                }
        }
    }


    void map (){
        // std::cout << "ENTER MAP\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        occurrencies.resize(itemsets.size());
        // hashmap.clear();
        //for every itemset
        size_t occ = 0;
        for (auto set : itemsets){
            occurrencies[occ]=0;
            //for every transaction
            for (auto& tx : transactions){
                unsigned int found = 0, cont=1;
                auto item = set.begin();
                //for every item in itemset
                while (cont && item != set.end() ){
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
                if (found == set.size()){
                    occurrencies[occ]++;
                }
            }
            ++occ;
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Map time: " << elapsed_seconds.count() << "s\n";

        // std::cout << "EXIT MAP\n";
    }

    void prune (float support){
        std::cout << "ITEMSETS BEFORE PRUNING: " << itemsets.size() << "\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();        
        size_t occ = 0;
        size_t size = transactions.size();
        auto item_set = std::begin(itemsets);
        while (item_set != std::end(itemsets)){
            if ((static_cast<float>(occurrencies[occ]) / (static_cast<float>(size))) < support){
                item_set = itemsets.erase(item_set);
            }
            else
                ++item_set;
            ++occ;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
        std::cout << "ITEMSETS AFTER PRUNING: " << itemsets.size() << "\n";
    }

    void merge (int k){
        // std::cout << "ENTER MERGE\n";

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        if (!itemsets.empty()){
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<std::vector<size_t>> temp;
            //for every itemset, try to unite it with another in the itemsets vector
            auto itemset_x = itemsets.begin();
            
            while (itemset_x!=itemsets.end()){
                auto itemset_y=itemset_x;
                itemset_y++;
                    while (itemset_y!=itemsets.end()){
                            std::vector<size_t> merged(k*k);
                            std::vector<size_t>::iterator it = std::set_union((*itemset_x).begin(), (*itemset_x).end(), (*itemset_y).begin(), (*itemset_y).end(), merged.begin());
                            merged.resize(it-merged.begin());
                            if (merged.size() == k)
                                temp.insert(merged);
                        itemset_y++;
                    }
                ++itemset_x;
            }
            itemsets = temp;
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout <<  "Merge time: " << elapsed_seconds.count() << "s\n";
    }


    void run (const std::string input_file, float support){
        int k=3;
        read_data(input_file);
        singles_merge(support);
        while (!itemsets.empty()){     
            std::cout << "ENTER PASS N° " << k << "\n";   
            map();
            prune(support);
            merge(k);
            std::cout << "EXIT PASS N° " << k << "\n\n";
            k++;
            
        }
    }
};