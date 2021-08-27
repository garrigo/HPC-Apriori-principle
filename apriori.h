#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_set>
#include <omp.h>

constexpr unsigned int byte_size = sizeof(unsigned int);
static constexpr size_t CACHE_SIZEE = 4194304;

struct VectorHash
{
    inline unsigned int operator()(const std::vector<unsigned int> &v) const
    {
        std::hash<unsigned int> hasher;
        unsigned int seed = 0;
        for (unsigned int i : v)
        {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

using VectorSet = std::unordered_set<std::vector<unsigned int>, VectorHash>;

struct Vector {
    typedef std::vector<std::vector<unsigned int>> data_type;
};

struct Set {
    typedef std::unordered_set<std::vector<unsigned int>, VectorHash> data_type;
};


template<class ItemsetType>
class AprioriBase
{
protected:
    std::vector<std::vector<unsigned int>> transactions;
    typename ItemsetType::data_type itemsets;
    std::vector<std::string> single_items;
    std::vector<unsigned int> occurrencies;
    size_t TX_BYTE_SIZE = 0;

    virtual void singles_merge(double, bool) = 0;
    virtual void map(unsigned int, bool) = 0;
    virtual void merge(unsigned int, double, bool) = 0;


    void store_itemsets(const std::string& filename)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if (ofs.is_open())
        {        
            for (auto& set : itemsets)
            {
                ofs << "{ ";
                for(auto& item: set)
                {
                    ofs << single_items[item] << " ";
                }
                ofs << "}";
                ofs << "  ";
            }
            ofs << "\n\n";
            ofs.close();
        }
        else std::cout << "Unable to open sse_output.dat file\n";
    }

    void read_data(const std::string input_file, bool parallel)
    {
        std::ifstream ifs;
        ifs.open(input_file);
        if (!ifs)
        {
            std::cout << "Input file could not be opened\n";
            exit(0);
        }
        std::string doc_buffer;

        #pragma omp parallel if(parallel)
        #pragma omp single
        while (!getline(ifs, doc_buffer).eof())
        {
            std::vector<unsigned int> line_buffer;
            std::istringstream iss(doc_buffer);
            std::string token;

            while (std::getline(iss, token, ' '))
            {
                if (!isspace(token[0]))
                {
                    bool clear = 1;
                    for (unsigned int i = 0; i < single_items.size(); i++)
                    {
                        if (single_items[i] == token)
                        {
                            clear = 0;
                            line_buffer.push_back(i);
                            occurrencies[i]++;
                            break;
                        }
                    }
                    if (clear)
                    {
                        single_items.push_back(token);
                        occurrencies.push_back(1);
                        line_buffer.push_back(single_items.size() - 1);
                    }
                }
            }
            #pragma omp task firstprivate(line_buffer), shared(transactions) if(parallel)
            {
                TX_BYTE_SIZE += line_buffer.size()*4;
                std::sort(line_buffer.begin(), line_buffer.end());
                #pragma omp critical (push)
                transactions.push_back(line_buffer);
            }
        }
        #pragma omp taskwait
        ifs.close();
    }

public:
    void run(const std::string input_file, double support, bool parallel)
    {
        int k = 2;
        read_data(input_file, parallel);
        singles_merge(support, parallel);
        while (!itemsets.empty())
        {
            map(k, parallel);
            ++k;
            merge(k, support, parallel);
        }
    }

};

class Apriori : public AprioriBase <Set>
{
    void singles_merge(double support, bool parallel)
    {
        #pragma omp parallel for schedule(dynamic) if(parallel)
        for (unsigned int i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned int j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    #pragma omp critical(singles_write)
                    {
                        itemsets.insert({i,j});
                    }
                }
        }
        // #pragma omp parallel if(parallel)
        // #pragma omp single
        // #pragma omp task
        // store_itemsets("nosse_set_output.dat");
    }

    void map(unsigned int k, bool parallel)
    {
        occurrencies.resize(itemsets.size());
        //for every itemset
        unsigned int occ = 0;
        #pragma omp parallel if(parallel)
        #pragma omp single
        for (const auto& set : itemsets)
        {        
            #pragma omp task firstprivate(set, occ) if(parallel) 
            {
                occurrencies[occ] = 0;
                //for every transaction
                for (const auto &tx : transactions)
                {
                    if (tx.size() >= k)
                    {
                        unsigned int found = 0, cont = 1, tx_cursor = 0;
                        auto item = set.begin();
                        //for every item in itemset
                        while (cont && item != set.end())
                        {
                            cont = 0;
                            //for every item in transaction
                            while (!cont && tx_cursor < tx.size())
                            {
                                if ((*item) < (tx[tx_cursor]))
                                {
                                    tx_cursor = tx.size();
                                }
                                else if ((*item) == (tx[tx_cursor]))
                                {
                                    ++found;
                                    cont = 1;
                                }

                                ++tx_cursor;
                            }
                            ++item;
                        }
                        if (found == k)
                        {
                            occurrencies[occ]++;
                        }
                    }
                }
            }
            ++occ;
        }
        #pragma omp taskwait
    }

    void prune(double support, bool parallel)
    {   
        unsigned int occ = 0;
        unsigned int size = transactions.size();
        auto item_set = std::begin(itemsets);
        while (item_set != std::end(itemsets))
        {
            if ((static_cast<double>(occurrencies[occ]) / (static_cast<double>(size))) < support)
            {
                item_set = itemsets.erase(item_set);
            }
            else
                ++item_set;
            ++occ;
        }
        // #pragma omp parallel if(parallel)
        // #pragma omp single
        // #pragma omp task
        // store_itemsets("nosse_set_output.dat");
    }

    void merge(unsigned int k, double support, bool parallel)
    {
        prune(support, parallel);
        if (!itemsets.empty())
        {
            VectorSet temp;
            //for every itemset, try to unite it with another in the itemsets vector
            auto itemset_x = itemsets.begin();
            #pragma omp parallel if(parallel)
            #pragma omp single
            while (itemset_x != itemsets.end())
            {
                    
                #pragma omp task firstprivate(itemset_x), shared(temp) if(parallel)
                {
                    auto itemset_y = itemset_x;
                    itemset_y++;
                    while (itemset_y != itemsets.end())
                    {
                        std::vector<unsigned int> merged(k);
                        unsigned int m = 0, v1 = 0, v2 = 0, distance = 0;
                        while (distance < 3 && m < k)
                        {
                            if (v2 == k - 1 || (v1 != k - 1 && (*itemset_x)[v1] < (*itemset_y)[v2]))
                            {
                                merged[m] = (*itemset_x)[v1];
                                ++v1;
                                ++distance;
                            }
                            else if (v1 == k - 1 || (v2 != k - 1 && (*itemset_x)[v1] > (*itemset_y)[v2]))
                            {
                                merged[m] = (*itemset_y)[v2];
                                ++v2;
                                ++distance;
                            }
                            else
                            {
                                merged[m] = (*itemset_y)[v2];
                                ++v1;
                                ++v2;
                            }
                            ++m;
                        }
                        if (distance == 2)
                            #pragma omp critical (merge_write)
                            {temp.insert(merged);}
                        itemset_y++;
                    }
                }

                ++itemset_x;
                
            }
            #pragma omp taskwait
            itemsets.swap(temp);
            // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }
    }
};


class SyncApriori : public AprioriBase<Vector>
{

    void singles_merge(double support, bool parallel)
    {
        const unsigned int tx_size = transactions.size(), single_size = single_items.size();
        // const unsigned int cache_regulator = std::max((500000/byte_size), 1U);
        #pragma parallel for schedule(dynamic) if(parallel)
        for (unsigned int i = 0; i < single_size - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(tx_size))) >= support)
                // #pragma omp for schedule(static) nowait
                for (unsigned int j = i + 1; j < single_size; j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(tx_size))) >= support)
                    {
                        itemsets.push_back({i, j});
                    }
                }
            // if (i%cache_regulator==0){
            //     #pragma omp barrier
            // }
        }
    }

    void map(unsigned int k, bool parallel)
    {
        
        occurrencies = std::vector<unsigned int>(itemsets.size(), 0);

        if(TX_BYTE_SIZE>(4*k*itemsets.size()))
        {
            unsigned int cache_regulator = CACHE_SIZEE/(TX_BYTE_SIZE/transactions.size());
            #pragma omp parallel if(parallel)
            for (unsigned int tx = 0; tx < transactions.size(); tx++)
            {
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned int set = 0; set < itemsets.size(); set++)
                {
                    unsigned int found = 0, cont = 1, tx_cursor = 0;
                    auto item = itemsets[set].begin();
                    //for every item in itemset
                    while (cont && item != itemsets[set].end())
                    {
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            if ((*item) < (transactions[tx][tx_cursor]))
                            {
                                break;
                            }
                            else if ((*item) == (transactions[tx][tx_cursor]))
                            {
                                ++found;
                                cont = 1;
                            }

                            ++tx_cursor;
                        }
                        ++item;
                    }
                    if (found == itemsets[set].size())
                    {
                        #pragma omp atomic
                        occurrencies[set]++;
                    }
                }
                if (tx % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
        else
        {
            unsigned int cache_regulator = CACHE_SIZEE/(k*4);
            #pragma omp parallel if(parallel)
            for (unsigned int set = 0; set < itemsets.size(); set++)
            {
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned int tx = 0; tx < transactions.size(); tx++)
                {
                    unsigned int found = 0, cont = 1, tx_cursor = 0;
                    auto item = itemsets[set].begin();
                    //for every item in itemset
                    while (cont && item != itemsets[set].end())
                    {
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            if ((*item) < (transactions[tx][tx_cursor]))
                            {
                                break;
                            }
                            else if ((*item) == (transactions[tx][tx_cursor]))
                            {
                                ++found;
                                cont = 1;
                            }

                            ++tx_cursor;
                        }
                        ++item;
                    }
                    if (found == itemsets[set].size())
                    {
                        #pragma omp atomic
                        occurrencies[set]++;
                    }
                }
                if (set % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
    }

    void merge(unsigned int k, double support, bool parallel)
    {

        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            VectorSet temp;
            std::vector<std::vector<unsigned int>> v_temp;
            unsigned int size = transactions.size();
            unsigned int itemsets_size = itemsets.size();
            size_t cache_regulator = CACHE_SIZEE/(k-1)*4;
            // const unsigned int cache_regulator = std::max((500000/byte_size*(k-1)), (unsigned int)1);
            //for every itemset, try to unite it with another in the itemsets vector

            // #pragma omp parallel
            #pragma omp parallel for schedule(dynamic) if(parallel)
            for (unsigned int i = 0; i < itemsets_size - 1; i++)
            {
                if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                {
                    // #pragma omp for schedule(static) nowait
                    for (unsigned int j = i + 1; j < itemsets_size; j++)
                    {
                        if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                        {
                            std::vector<unsigned int> merged(k);
                            unsigned int m = 0, v1 = 0, v2 = 0, distance = 0;
                            while (distance < 3 && m < k)
                            {
                                if (v2 == k - 1 || (v1 != k - 1 && itemsets[i][v1] < itemsets[j][v2]))
                                {
                                    merged[m] = itemsets[i][v1];
                                    ++v1;
                                    ++distance;
                                }
                                else if (v1 == k - 1 || (v2 != k - 1 && itemsets[i][v1] > itemsets[j][v2]))
                                {
                                    merged[m] = itemsets[j][v2];
                                    ++v2;
                                    ++distance;
                                }
                                else
                                {
                                    merged[m] = itemsets[j][v2];
                                    ++v1;
                                    ++v2;
                                }
                                ++m;
                            }
                            if (distance == 2)
                            {
                                #pragma omp critical(merge_write)
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
            // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }
    }
};