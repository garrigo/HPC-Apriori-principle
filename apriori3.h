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

class Apriori
{
public:
    std::vector<std::vector<std::string>> transactions;
    std::vector<std::set<std::string>> itemsets;
    std::unordered_set<std::string> single_items;
    std::unordered_map<std::string, size_t> hashmap;

    inline void print_single_items()
    {
        for (auto i : single_items)
            std::cout << "Element: " << i << " - Cardinality: " << hashmap[i] << "\n";
        std::cout << "SINGLE ITEMS Total size: " << single_items.size() << "\n";
    }

    inline void print_items()
    {
        for (auto set : itemsets)
        {
            for (auto str : set)
            {
                std::cout << str << "    ";
            }
            std::cout << "\n";
        }
        std::cout << "ITEMSETS Total size: " << itemsets.size() << "\n";
    }

    inline void print_items(std::set<std::string> &set, std::string &buff)
    {
        for (auto str : set)
        {
            std::cout << str << "    ";
        }
        std::cout << " - Cardinality: " << static_cast<double>(hashmap[buff]) / static_cast<double>(transactions.size()) << "\n";
    }

    inline void update_single_structures(std::string token)
    {
        single_items.insert(token);
        hashmap[token]++;
    }

    void read_data(const std::string input_file)
    {
        std::ifstream ifs;
        ifs.open(input_file);
        if (!ifs)
        {
            std::cout << "Input file could not be opened\n";
            exit(0);
        }
        std::string doc_buffer;
        std::vector<std::vector<std::string>> result;

        while (!getline(ifs, doc_buffer).eof())
        {
            std::vector<std::string> line_buffer;
            std::istringstream iss(doc_buffer);
            std::string token;
            while (std::getline(iss, token, ' '))
            {
                if (!isspace(token[0]))
                {
                    line_buffer.push_back(token);
                    update_single_structures(token);
                }
            }
            transactions.push_back(line_buffer);
        }
    }

    void singles_prune(double support)
    {
        auto item = std::begin(single_items);
        while (item != std::end(single_items))
        {
            // std::cout << *item << " Cardinality: " << hashmap[*item]/transactions.size() << "\n";
            if (static_cast<double>(hashmap[*item]) / static_cast<double>(transactions.size()) < support)
                item = single_items.erase(item);
            else
                ++item;
        }
    }

    void singles_merge()
    {
        auto item_x = std::begin(single_items);
        auto item_y = item_x;
        while (item_x != std::end(single_items))
        {
            ++item_y;
            while (item_y != std::end(single_items))
            {
                itemsets.push_back({*item_x, *item_y});
                item_y++;
            }
            ++item_x;
            item_y = item_x;
        }
    }

    // void singles_merge(){
    //     auto item_x = std::begin (single_items);
    //     auto item_y = item_x;
    //     while (item_x != std::end(single_items)){
    //         ++item_y;
    //         while (item_y != std::end(single_items)){
    //             itemsets.push_back({*item_x, *item_y});
    //             item_y++;
    //         }
    //         ++item_x;
    //         item_y = item_x;
    //     }
    // }

    void map()
    {
        // std::cout << "ENTER MAP\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        hashmap.clear();
        //for every itemset
        for (auto &itemset : itemsets)
        {
            std::string buff;
            //for every transaction
            for (auto &tx : transactions)
            {
                buff.clear();
                unsigned int found = 0, cont = 1;
                auto item = std::begin(itemset);
                //for every item in itemset
                while (cont && item != std::end(itemset))
                {
                    buff += (*item);
                    cont = 0;
                    //for every item in transaction
                    for (auto &item_tx : tx)
                    {
                        if (*item == item_tx)
                        {
                            found++;
                            cont = 1;
                            break;
                        }
                    }
                    item++;
                }
                if (found == itemset.size())
                {
                    hashmap[buff]++;
                }
            }
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        // std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
        // std::cout << "EXIT MAP\n";
    }

    /*     void map (){
        // std::cout << "ENTER MAP\n";
        hashmap.clear();
        //for every transaction
        for (auto& tx : transactions){
            //for every itemset
            for (auto& itemset : itemsets){
                std::string buff;
                unsigned int found = 0, cont=1;
                auto item = std::begin (itemset);
                //for every item in itemset
                while (cont && item != std::end(itemset) ){
                    //for every item in transaction
                    buff+=(*item);
                    cont = 0;
                    for (auto& item_tx : tx){
                        if (*item==item_tx){
                            found++;
                            cont = 1;
                            break;
                        }
                    }
                    item++;
                }
                if (found == itemset.size()){
                    hashmap[buff]++;
                }

            }
        }
        // std::cout << "EXIT MAP\n";
    } */

    void prune(double support)
    {
        //PRUNE
        // std::cout << "ENTER PRUNE\n";
        // std::cout << "ITEMSETS BEFORE PRUNING: " << itemsets.size() << "\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        auto itemset = std::begin(itemsets);
        while (itemset != std::end(itemsets))
        {
            std::string buff;
            auto item = std::begin(*itemset);
            while (item != std::end(*itemset))
            {
                buff += (*item);
                item++;
            }
            // print_items(*itemset, buff);
            if (static_cast<double>(hashmap[buff]) / static_cast<double>(transactions.size()) < support)
                itemset = itemsets.erase(itemset);
            else
                ++itemset;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        // std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
        // std::cout << "ITEMSETS AFTER PRUNING: " << itemsets.size() << "\n";
        // std::cout << "EXIT PRUNE\n";
    }

    //     void merge (int k){
    //         std::cout << "ENTER MERGE\n";
    //             std::chrono::time_point<std::chrono::system_clock> start, end;
    //     start = std::chrono::system_clock::now();
    //         if (!itemsets.empty()){
    //             // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
    //             std::set<std::set<std::string>> temp;
    //             //for every itemset, try to unite it with another in the itemsets vector
    //             for (int itemset_x=0; itemset_x<itemsets.size(); itemset_x++){
    //                 auto cursor = single_items.begin();
    //                 for (int single=0; single<single_items.size(); single++){
    //                     std::set<std::string> merged(itemsets[itemset_x]);
    //                     merged.insert(*cursor);
    //                     cursor++;
    //                     if (merged.size() == k)
    //                         temp.insert(merged);
    //                 }
    //             }
    //             itemsets.assign( temp.begin(), temp.end() );
    //         }
    //         end = std::chrono::system_clock::now();
    //         std::chrono::duration<double> elapsed_seconds = end - start;
    // std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
    //         std::cout << "EXIT MERGE PASS N° " << k << "\n";
    //     }

    void merge(int k)
    {
        // std::cout << "ENTER MERGE\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<std::set<std::string>> temp;
            //for every itemset, try to unite it with another in the itemsets vector
            for (int itemset_x = 0; itemset_x < itemsets.size() - 1; itemset_x++)
            {
                for (int itemset_y = itemset_x + 1; itemset_y < itemsets.size(); itemset_y++)
                {
                    std::set<std::string> merged(itemsets[itemset_x]);
                    merged.insert(itemsets[itemset_y].begin(), itemsets[itemset_y].end());
                    if (merged.size() == k)
                        temp.insert(merged);
                }
            }
            itemsets.assign(temp.begin(), temp.end());
        }
        std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        // std::cout <<  "Elapsed time: " << elapsed_seconds.count() << "s\n";
    }

    void run(const std::string input_file, double support)
    {
        int k = 3;
        read_data(input_file);
        singles_prune(support);
        singles_merge();
        // print_single_items();
        // print_items();
        while (!itemsets.empty())
        {
            // std::cout << "ENTER PASS N° " << k << "\n\n";
            map();
            prune(support);
            merge(k);
            // std::cout << "EXIT PASS N° " << k << "\n\n";
            k++;
        }
    }
};
