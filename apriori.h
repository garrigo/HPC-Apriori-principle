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



struct VectorHash
{
    inline unsigned operator()(const std::vector<unsigned> &v) const
    {
        std::hash<unsigned> hasher;
        unsigned seed = 0;
        for (unsigned i : v)
        {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

using VectorSet = std::unordered_set<std::vector<unsigned>, VectorHash>;

struct Vector {
    typedef std::vector<unsigned>data_type;
};

struct Set {
    typedef std::unordered_set<std::vector<unsigned>, VectorHash> data_type;
};


template<class ItemsetType>
class AprioriBase
{
protected:
    std::vector<std::vector<unsigned>> transactions;
    typename ItemsetType::data_type itemsets;
    std::vector<std::string> single_items;
    std::vector<unsigned> occurrencies;
    size_t TX_BYTE_SIZE = 0;
    static constexpr size_t CACHE_SIZE = 1048576;

    virtual void singles_merge(const double, const int) = 0;
    virtual void map(const unsigned, const int) = 0;
    virtual void merge(const unsigned, const double, const int) = 0;

    //Output single items vector to file
    void store_single_items(std::string filename, double support)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if(ofs.is_open())
        {
            for(unsigned i = 0; i<single_items.size(); i++)
            {
                if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                    ofs << "{" << single_items[i] << "} ";
            }
            ofs << "\n\n";
            ofs.close();
        }
        else std::cout << "Unable to open sse_output.dat file\n";

    }

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

    void read_data(const std::string input_file, const int max_threads)
    {
        std::ifstream ifs;
        ifs.open(input_file);
        if (!ifs)
        {
            std::cout << "Input file could not be opened\n";
            exit(0);
        }
        std::string doc_buffer;

        #pragma omp parallel num_threads(max_threads)
        #pragma omp single
        while (!getline(ifs, doc_buffer).eof())
        #pragma omp task firstprivate(doc_buffer)
        {
            std::vector<unsigned> line_buffer;
            std::istringstream iss(doc_buffer);
            std::string token;

            while (std::getline(iss, token, ' '))
            {
                if (!isspace(token[0]))
                {
                    bool clear = 1;
                    unsigned i, temp_size;
                    #pragma omp critical (push_single)
                    temp_size = single_items.size();
                    // #pragma omp critical (push_single)
                    {
                        for (i = 0; i < temp_size; i++)
                        {
                            if (single_items[i] == token)
                            {
                                clear = 0;
                                #pragma omp atomic
                                occurrencies[i]++;
                                break;
                            }
                        }
                        if (clear)
                        {
                            #pragma omp critical (push_single)
                            {
                                if (temp_size != single_items.size())
                                {
                                    for (i = temp_size; i < single_items.size(); i++)
                                    {
                                        if (single_items[i] == token)
                                        {
                                            clear = 0;
                                            #pragma omp atomic
                                            occurrencies[i]++;
                                            break;
                                        }
                                    }                                    
                                }
                                if(clear){
                                    single_items.push_back(token);
                                    occurrencies.push_back(1); 
                                }                               
                            }
                        }
                    }
                    line_buffer.push_back(i);
                }
            }
            #pragma omp atomic
            TX_BYTE_SIZE += line_buffer.size()*4;
            std::sort(line_buffer.begin(), line_buffer.end());
            #pragma omp critical (push_transactions)
            transactions.push_back(line_buffer);
            
        }
        ifs.close();
    }

public:
    void run(const std::string input_file, const double support, const int max_threads)
    {
        unsigned k = 2;
        read_data(input_file, max_threads);
        singles_merge(support, max_threads);
        while (!itemsets.empty())
        {
            map(k, max_threads);
            ++(k);
            merge(k, support, max_threads);
        }
    }

};

class Apriori : public AprioriBase <Set>
{
    void singles_merge(const double support, const int max_threads)
    {
        #pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads)
        for (unsigned i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    #pragma omp critical(singles_write)
                    {
                        itemsets.insert({i,j});
                    }
                }
        }
        // store_single_items("apriori_out.dat", support);
        // store_itemsets("nosse_set_output.dat");
    }

    void map(const unsigned k, const int max_threads)
    {
        occurrencies.resize(itemsets.size());
        //for every itemset
        unsigned occ = 0;
        #pragma omp parallel num_threads(max_threads)
        #pragma omp single
        for (const auto& set : itemsets)
        {        
            #pragma omp task firstprivate(set, occ)
            {
                occurrencies[occ] = 0;
                //for every transaction
                for (const auto &tx : transactions)
                {
                    if (tx.size() >= k)
                    {
                        unsigned found = 0, cont = 1, tx_cursor = 0;
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

    void prune(const double support, const int max_threads)
    {   
        unsigned occ = 0;
        unsigned size = transactions.size();
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
        // store_itemsets("nosse_set_output.dat");
    }

    void merge(const unsigned k, const double support, const int max_threads)
    {
        prune(support, max_threads);
        if (!itemsets.empty())
        {
            VectorSet temp;
            //for every itemset, try to unite it with another in the itemsets vector
            auto itemset_x = itemsets.begin();
            #pragma omp parallel num_threads(max_threads)
            #pragma omp single
            while (itemset_x != itemsets.end())
            {
                    
                #pragma omp task firstprivate(itemset_x), shared(temp)
                {
                    auto itemset_y = itemset_x;
                    itemset_y++;
                    while (itemset_y != itemsets.end())
                    {
                        std::vector<unsigned> merged(k);
                        unsigned m = 0, v1 = 0, v2 = 0, distance = 0;
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
            itemsets.swap(temp);
            // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }
    }
};


class SyncApriori : public AprioriBase<Vector>
{

    void store_itemsets(const std::string& filename, unsigned k)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if (ofs.is_open())
        {      
            for (unsigned i = 0; i<itemsets.size(); i+=k)
            {
                ofs << "{ ";
                for(unsigned j = 0; j<k; j++)
                {
                    ofs << single_items[itemsets[i+j]] << " ";
                }
                ofs << "}";
                ofs << "  ";
            }
            ofs << "\n\n";
            ofs.close();
        }
        else std::cout << "Unable to open sse_output.dat file\n";
    }

    void singles_merge(const double support, const int max_threads)
    {
        const unsigned cache_regulator = CACHE_SIZE/8;
        const unsigned tx_size = transactions.size(), 
                            single_size = single_items.size();
        #pragma omp parallel num_threads(max_threads)
        #pragma omp for schedule(dynamic,1)
        for (unsigned i = 0; i < single_size - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(tx_size))) >= support)
                // #pragma omp for schedule(static) nowait 
                for (unsigned j = i + 1; j < single_size; j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(tx_size))) >= support)
                    #pragma omp critical (single_write)
                    {
                        itemsets.push_back(i);
                        itemsets.push_back(j);
                    }
                }
            // if (i % cache_regulator == 0)
            // {
            //     #pragma omp barrier
            // } 
        }
        // store_single_items("sync_out.dat", support);
        // store_itemsets("sync_out.dat", 2);

    }

    void map(const unsigned k, const int max_threads)
    {
        const unsigned itemsets_size = itemsets.size()/k;
        occurrencies = std::vector<unsigned>(itemsets_size, 0);
        

        if(TX_BYTE_SIZE>(4*itemsets.size()))
        {
            const unsigned cache_regulator = CACHE_SIZE/(TX_BYTE_SIZE/transactions.size());
            #pragma omp parallel num_threads(max_threads)
            for (unsigned tx = 0; tx < transactions.size(); tx++)
            {
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned set = 0; set < itemsets_size; set++)
                {
                    unsigned found = 0, cont = 1, tx_cursor = 0, item = 0;
                    //for every item in itemset
                    while (cont && item < k)
                    {
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            if (itemsets[set*k+item] < (transactions[tx][tx_cursor]))
                            {
                                break;
                            }
                            else if (itemsets[set*k+item] == (transactions[tx][tx_cursor]))
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
            const unsigned cache_regulator = CACHE_SIZE/(k*4);
            #pragma omp parallel num_threads(max_threads)
            for (unsigned set = 0; set < itemsets_size; set++)
            {
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned tx = 0; tx < transactions.size(); tx++)
                {
                    unsigned found = 0, cont = 1, tx_cursor = 0, item = 0;
                    //for every item in itemset
                    while (cont && item < k)
                    {
                        cont = 0;
                        //for every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            if (itemsets[set*k+item] < (transactions[tx][tx_cursor]))
                            {
                                break;
                            }
                            else if (itemsets[set*k+item] == (transactions[tx][tx_cursor]))
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

    void prune(const double support, const unsigned k, const int max_threads)
    {   
        unsigned size = transactions.size();
        // unsigned tail = itemsets.size() - 1;
        std::vector<unsigned> temp;
        #pragma omp parallel num_threads(max_threads)
        {
            #pragma omp for schedule(static) 
            for(int i = 0; i<itemsets.size(); i+=k)
            {
                if ((static_cast<double>(occurrencies[i/k]) / (static_cast<double>(size))) >= support)
                {   
                    
                    #pragma omp critical(prune)
                    temp.insert(temp.end(), itemsets.begin()+i, itemsets.begin()+i+k);
                }
            }

        }
        itemsets.swap(temp);
        // store_itemsets("sync_out.dat", k);
    }


    void merge(const unsigned k, const  double support, const int max_threads)
    {
        prune(support, k-1, max_threads);
        if (!itemsets.empty())
        {
            std::set<std::vector<unsigned>> temp;
            const unsigned size = transactions.size();
            const unsigned itemsets_size = itemsets.size()/(k-1);
            const unsigned cache_regulator = CACHE_SIZE/(k*4);

            #pragma omp parallel num_threads(max_threads)
            {
                #pragma omp for schedule(dynamic, 1)
                for (unsigned i = 0; i < itemsets_size - 1; i++)
                {
                    // if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                    {
                        unsigned i_off = i*(k-1);
                        // #pragma omp for schedule(static) nowait 
                        for (unsigned j = i + 1; j < itemsets_size; j++)
                        {

                            // if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                            {
                                
                                std::vector<unsigned> merged(k);
                                unsigned m = 0, v1 = 0, v2 = 0, distance = 0, j_off = j*(k-1);
                                while (distance < 3 && m < k)
                                {
                                    if (v2 == k - 1 || (v1 != k - 1 && itemsets[i_off+v1] < itemsets[j_off+v2]))
                                    {
                                        merged[m] = itemsets[i_off+v1];
                                        ++v1;
                                        ++distance;
                                    }
                                    else if (v1 == k - 1 || (v2 != k - 1 && itemsets[i_off+v1] > itemsets[j_off+v2]))
                                    {
                                        merged[m] = itemsets[j_off+v2];
                                        ++v2;
                                        ++distance;
                                    }
                                    else
                                    {
                                        merged[m] = itemsets[j_off+v2];
                                        ++v1;
                                        ++v2;
                                    }
                                    ++m;
                                }
                                if (distance == 2)
                                {
                                    #pragma omp critical(merge_write)
                                    temp.insert(merged);
                                    
                                }
                            }
                        }
                    }
                    // if (i % cache_regulator == 0)
                    // {
                    //     #pragma omp barrier
                    // }
                }
                // #pragma omp barrier
                #pragma omp single
                {
                    itemsets.resize(temp.size()*k);
                    unsigned i = 0;
                    for (auto& set : temp)
                    {
                        #pragma omp task firstprivate(set, i), shared(itemsets)
                        for (unsigned j = 0; j<k; j++)
                            itemsets[i+j] = set[j];
                        i+=k;
                        
                    }
                } 
            }
            // std::cout << "ITEMSETS SIZE: " << itemsets.size()/k << "\n";
        }
    }
};