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


//Hash operation to allocate a vector in a std::unordered_set
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

//Data type of the itemset collection for the Sync version
struct Vector {
    typedef std::vector<unsigned>data_type;
};

//Data type of the itemset collection for the normal version
struct Set {
    typedef std::unordered_set<std::vector<unsigned>, VectorHash> data_type;
};


//Abstract class that includes all members that are shared by all No-SSE versions
template<class ItemsetType>
class AprioriBase
{
protected:
    //Vector of transactions taken from the file
    std::vector<std::vector<unsigned>> transactions;
    //Collection of itemsets generated at every step
    typename ItemsetType::data_type itemsets;
    //Vector of unique items taken by the file
    std::vector<std::string> single_items;
    //Vector representing the occurrencies of every itemset in the whole transaction list
    std::vector<unsigned> occurrencies;

    //Counter for the average number of bytes of a transaction line
    size_t TX_BYTE_SIZE = 0;
    

    //Methods to be implemented by every child class

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

    //Output current itemsets to file
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

    //Read the input file and populate single_items, transactions and occurrencies
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
        //Create a task for every line of the input file
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
                    //Search for the new token in the vector of single items already found (no-critical section pass)
                    for (i = 0; i < temp_size; i++)
                    {
                        if (single_items[i] == token)
                        {
                            clear = 0;
                            //If token is found in the vector, increment atomically the corresponding entry of occurrencies
                            #pragma omp atomic
                            occurrencies[i]++;
                            break;
                        }
                    }
                    //Token wasn't in the vector during the no-critical section search
                    if (clear)
                    {
                        /*Try again to find the token in another pass,
                        but this time in a critical section and by starting from the last element compared in the first search*/                        
                        #pragma omp critical (push_single)
                        {
                            if (temp_size != single_items.size())
                            {
                                for (i = temp_size; i < single_items.size(); i++)
                                {
                                    if (single_items[i] == token)
                                    {
                                        clear = 0;
                                        //If token is found in the vector, increment atomically the corresponding entry of occurrencies
                                        #pragma omp atomic
                                        occurrencies[i]++;
                                        break;
                                    }
                                }                                    
                            }
                            //If the token wasn't found again, create a new entry both in the single items vector and in the occurrencies' one
                            if(clear){
                                single_items.push_back(token);
                                occurrencies.push_back(1); 
                            }                               
                        }
                    }
                    //Save the position of the token (in the single_items vector) in the new transaction line
                    line_buffer.push_back(i);
                }
            }
            #pragma omp atomic
            TX_BYTE_SIZE += line_buffer.size()*4;
            //Sort the new transaction line and push it in the vector of transactions
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
        //While itemsets within the given threshold are being found
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
    //Create itemsets by merging couples of single items (for the k=2 pass)
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

    //For every itemset, search it in every transaction and if it's present increment the corresponding occurrencies' entry
    void map(const unsigned k, const int max_threads)
    {
        occurrencies.resize(itemsets.size());
        unsigned occ = 0;
        #pragma omp parallel num_threads(max_threads)
        #pragma omp single
        //For every itemset create a task
        for (const auto& set : itemsets)
        {        
            #pragma omp task firstprivate(set, occ)
            {
                occurrencies[occ] = 0;
                //For every transaction
                for (const auto &tx : transactions)
                {
                    if (tx.size() >= k)
                    {
                        bool cont = true;
                        unsigned found = 0, tx_cursor = 0;
                        auto item = set.begin();
                        //For every item in itemset
                        while (cont && item != set.end())
                        {
                            cont = false;
                            //For every item in transaction
                            while (!cont && tx_cursor < tx.size())
                            {
                                /*If the item of the itemset is less than the item of the transaction, it means that the item is not present in that transaction line,
                                since both the itemset and the transaction are sorted in ascending order. In this case we can skip to the next transaction.*/
                                if ((*item) < (tx[tx_cursor]))
                                    break;
                                else if ((*item) == (tx[tx_cursor]))
                                {
                                    ++found;
                                    cont = true;
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

    //Prune all itemsets that do not reach the threshold
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

    //Merge all k-1 pass itemsets to create a new set of itemsets having k single items
    void merge(const unsigned k, const double support, const int max_threads)
    {
        prune(support, max_threads);
        if (!itemsets.empty())
        {
            VectorSet temp;
            
            auto itemset_x = itemsets.begin();
            #pragma omp parallel num_threads(max_threads)
            #pragma omp single
            //A single thread creates a task for every itemset
            while (itemset_x != itemsets.end())
            {
                //The task tries to unite the given itemset with all the itemsets following it
                #pragma omp task firstprivate(itemset_x), shared(temp)
                {
                    auto itemset_y = itemset_x;
                    itemset_y++;
                    while (itemset_y != itemsets.end())
                    {
                        std::vector<unsigned> merged(k);
                        unsigned m = 0, v1 = 0, v2 = 0, distance = 0;
                        /*Since items in itemsets are sorted, compare one by one (starting from the beginning) the elements of the two itemsets.
                        Two itemsets of size (k-1) merge into one of size k only if they differ just by one element. By comparing one by one,
                        we can stop the merging when the "distance" between the two vectors reaches value 3, which means that there are too many
                        different items to generate a k-vector.*/
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
                        //If the merged set has size k try to put it in the new set of itemsets
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
    //Custom cache-aware value to compute a cache_regulator
    static constexpr size_t CACHE_SIZE = 1048576;

    //Output current itemsets to file
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

    //Create itemsets by merging couples of single items (for the k=2 pass)
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

    //For every itemset, search it in every transaction and if it's present increment the corresponding occurrencies' entry
    void map(const unsigned k, const int max_threads)
    {
        const unsigned itemsets_size = itemsets.size()/k;
        occurrencies = std::vector<unsigned>(itemsets_size, 0);

        if(TX_BYTE_SIZE>(4*itemsets.size()))
        {
            const unsigned cache_regulator = CACHE_SIZE/(TX_BYTE_SIZE/transactions.size());
            #pragma omp parallel num_threads(max_threads)
            //For every transaction
            for (unsigned tx = 0; tx < transactions.size(); tx++)
            {
                //For every itemset
                #pragma omp for schedule(static) nowait
                for (unsigned set = 0; set < itemsets_size; set++)
                {
                    bool cont = true;
                    unsigned found = 0, tx_cursor = 0, item = 0;
                    //For every item in itemset
                    while (cont && item < k)
                    {
                        cont = false;
                        //For every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            
                            /*If the item of the itemset is less than the item of the transaction, it means that the item is not present in that transaction line,
                            since both the itemset and the transaction are sorted in ascending order. In this case we can skip to the next transaction.*/
                            if (itemsets[set*k+item] < (transactions[tx][tx_cursor]))
                                break;
                            else if (itemsets[set*k+item] == (transactions[tx][tx_cursor]))
                            {
                                ++found;
                                cont = true;
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
                    bool cont = true;
                    unsigned found = 0, tx_cursor = 0, item = 0;
                    //for every item in itemset
                    while (cont && item < k)
                    {
                        cont = false;
                        //for every item in transaction
                        while (!cont && tx_cursor < transactions[tx].size())
                        {
                            /*If the item of the itemset is less than the item of the transaction, it means that the item is not present in that transaction line,
                            since both the itemset and the transaction are sorted in ascending order. In this case we can skip to the next transaction.*/
                            if (itemsets[set*k+item] < (transactions[tx][tx_cursor]))
                                break;
                            else if (itemsets[set*k+item] == (transactions[tx][tx_cursor]))
                            {
                                ++found;
                                cont = true;
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

    //Prune all itemsets that do not reach the threshold
    void prune(const double support, const unsigned k, const int max_threads)
    {   
        unsigned size = transactions.size();
        std::vector<unsigned> temp;
        #pragma omp parallel num_threads(max_threads)
        {
            #pragma omp for schedule(dynamic, 1) 
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

    //Merge all k-1 pass itemsets to create a new set of itemsets having k single items
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
                                /*Since items in itemsets are sorted, compare one by one (starting from the beginning) the elements of the two itemsets.
                                Two itemsets of size (k-1) merge into one of size k only if they differ just by one element. By comparing one by one,
                                we can stop the merging when the "distance" between the two vectors reaches value 3, which means that there are too many
                                different items to generate a k-vector.*/
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
                                //If the merged set has size k try to put it in the new set of itemsets
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
                    //Fill the itemsets vector with the new itemsets using the temporary set of the merge phase
                    //For every set create a task that fills the unidimensional vector
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