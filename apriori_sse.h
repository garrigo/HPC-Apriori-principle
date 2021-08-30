#pragma once
#include <vector>
#include <math.h>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <omp.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <bitset>

static unsigned int MASK_SIZE;
static constexpr size_t CACHE_SIZE = 4194304;

struct Compare
{
    bool operator()(unsigned int * __restrict__ a, unsigned int * __restrict__ b) const
    {
        unsigned int compare_result[8];
        __m128i *p_a = (__m128i *)a, *p_b = (__m128i *)b;
        for (unsigned int mask = 0; mask < MASK_SIZE / 128; mask++)
        {
            __m128i m_a = _mm_load_si128(p_a);
            __m128i m_b = _mm_load_si128(p_b);
            _mm_storeu_si128((__m128i *)compare_result, _mm_cmpgt_epi32(m_b, m_a));
            _mm_storeu_si128((__m128i *)(compare_result + 4), _mm_cmpgt_epi32(m_a, m_b));
            for (unsigned int b = 0; b < 4; b++)
            {
                if (compare_result[b])
                    return true;
                if (compare_result[b + 4])
                    return false;
            }
            ++p_a;
            ++p_b;
        }
        return false;
    }
};


struct VectorSSE {
    typedef std::vector<unsigned int *> data_type;
};

struct SetSSE {
    typedef std::set<unsigned int *, Compare> data_type;
};


template<class ItemsetType>
class SSE
{
protected:
    std::vector<std::string> single_items;
    std::vector<unsigned int *> transactions;
    std::vector<unsigned int> masks;
    typename ItemsetType::data_type  itemsets;
    std::vector<unsigned int> occurrencies;

    virtual void singles_merge(const double, const int) = 0;
    virtual void map(const unsigned int, const int) = 0;
    virtual void merge(const unsigned int, const double, const int) = 0;
  
    void store_single_items(std::string filename)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if(ofs.is_open())
        {
            for(auto& s : single_items)
            {
                ofs << s << " ";
            }
            ofs << "\n\n";
            ofs.close();
        }
        else std::cout << "Unable to open sse_output.dat file\n";

    }

    void store_itemsets(std::string filename)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if (ofs.is_open())
        {        
            for (auto& set : itemsets)
            {
                for (unsigned int i = 0; i < MASK_SIZE / 32; i++)
                {
                    std::bitset<32> x(set[i]);
                    ofs << x << " ";
                }
                ofs << "\n";
            }
            ofs << "\n";
            ofs.close();
        }
        else std::cout << "Unable to open sse_output.dat file\n";
    }

    inline void initialize_array(unsigned int *buf)
    {
        for (unsigned int m = 0; m < MASK_SIZE / 128; m++)
            _mm_store_si128((__m128i *)(buf + (m * 4)), _mm_setzero_si128());
    }

    void read_data(const std::string input_file, int max_threads)
    {

        std::ifstream ifs;
        ifs.open(input_file);
        if (!ifs)
        {
            std::cout << "Input file could not be opened\n";
            exit(0);
        }
        std::string doc_buffer;
        while (!getline(ifs, doc_buffer).eof())
        {
            std::istringstream iss(doc_buffer);
            std::string token;
            unsigned int *buf = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);
            initialize_array(buf);
            while (std::getline(iss, token, ' '))
            {

                if (!isspace(token[0]))
                {
                    bool clear = true;
                    for (unsigned int i = 0; i < single_items.size(); i++)
                    {
                        if (single_items[i] == token)
                        {
                            clear = 0;
                            occurrencies[i]++;
                            buf[i / 32] ^= static_cast<unsigned int>(std::exp2(31 - (i % 32)));
                            // i=single_items.size();
                            break;
                        }
                    }
                    //new element found
                    if (clear)
                    {
                        single_items.push_back(token);
                        occurrencies.push_back(1);
                        unsigned int s_size = single_items.size();
                        if (single_items.size() > MASK_SIZE)
                        {
                            MASK_SIZE += 128;
                            unsigned int *buf2 = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);

                            for (unsigned int i = 0; i < MASK_SIZE / 128 - 1; i++)
                                _mm_store_si128((__m128i *)(buf2 + (i * 4)), _mm_load_si128((__m128i *)(buf + (i * 4))));
                            _mm_store_si128((__m128i *)(buf2 + MASK_SIZE / 32 - 4), _mm_setzero_si128());
                            _mm_free(buf);
                            buf = buf2;
                        }

                        buf[(s_size - 1) / 32] ^= static_cast<unsigned int>(std::exp2(31 - ((s_size - 1) % 32) ));
                    }
                }
            }
            transactions.push_back(buf);
            masks.push_back(MASK_SIZE);
        }
        ifs.close();
    } 

public:
    void run(const std::string input_file, const double support, const int max_threads)
    {   
        
        unsigned int k = 2;
        MASK_SIZE = 128;
        read_data(input_file, max_threads);
        singles_merge(support, max_threads); 
        while (!itemsets.empty())
        {
            map(k, max_threads);
            ++k;
            merge(k, support, max_threads);
        }
        #pragma omp parallel for num_threads(max_threads)
        for (unsigned int i = 0; i < transactions.size(); i++)
            _mm_free(transactions[i]);
    }
};



class AprioriSSE : public SSE<SetSSE>
{

    void singles_merge(const double support, const int max_threads)
    {
        #pragma omp parallel for num_threads(max_threads) schedule(dynamic)
        for (unsigned int i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned int j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        unsigned int *buf = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);
                        initialize_array(buf);
                        buf[i / 32] ^= static_cast<unsigned int>(std::exp2(31 - (i % 32) ));
                        buf[j / 32] ^= static_cast<unsigned int>(std::exp2(31 - (j % 32) ));
                        #pragma omp critical(single_write)
                        {
                            itemsets.insert(buf);
                        }
                    }
                }
        }
        // #pragma omp parallel num_threads(max_threads)
        // #pragma single
        // #pragma omp task
        // {        
        //     store_single_items("sse_output.dat");
        //     store_itemsets("sse_output.dat");
        // }
    }

    void map(const unsigned int k, const int max_threads)
    {
        
        if(transactions.size() < itemsets.size())
        {
            unsigned int occ = 0;
            occurrencies.resize(itemsets.size());
            //for every itemset
            #pragma omp parallel num_threads(max_threads)
            #pragma omp single
            for (const auto set : itemsets)
            {
                occurrencies[occ] = 0;
                //for every transaction
                #pragma omp task firstprivate(set, occ)
                for (unsigned int tx = 0; tx < transactions.size(); tx++)
                {
                    bool found = true;
                    unsigned int masks_number = (masks[tx] / 128);
                    __m128i *p_set = (__m128i *)(set), *p_tx = (__m128i *)transactions[tx];
                    for (unsigned int i = 0; i < MASK_SIZE / 128; i++)
                    {

                        __m128i set_128 = _mm_load_si128(p_set);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(p_tx)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                        ++p_set;
                        ++p_tx;
                    }
                    if (found)
                        occurrencies[occ]++;
                }
                ++occ;
            }
            #pragma omp taskwait
        }
        else
        {
            occurrencies = std::vector<unsigned int>(itemsets.size(), 0);

            #pragma omp parallel for schedule(static) num_threads(max_threads)
            for (unsigned int tx = 0; tx < transactions.size(); tx++)
            {
                unsigned int occ = 0;
                unsigned int masks_number = (masks[tx] / 128);
                for (const auto set : itemsets)
                {
                    bool found = true;
                    __m128i *p_set = (__m128i *)(set), *p_tx = (__m128i *)transactions[tx];
                    for (unsigned int i = 0; i < MASK_SIZE / 128; i++)
                    {

                        __m128i set_128 = _mm_load_si128(p_set);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(p_tx)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                        ++p_set;
                        ++p_tx;
                    }
                    if (found)
                        #pragma omp atomic
                        occurrencies[occ]++;
                    ++occ;
                }
            }
        }
    }

    void prune(const double support, const int max_threads)
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
        // #pragma omp parallel num_threads(max_threads)
        // #pragma single
        // #pragma omp task
        // {        
        //     store_itemsets("sse_output.dat");
        // }
    }

    void merge(const unsigned int k, const double support, const int max_threads)
    {
        prune(support, max_threads);
        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<unsigned int *, Compare> temp;
            auto itemset_x = itemsets.begin();
            //for every itemset, try to unite it with another in the itemsets vector
            #pragma omp parallel num_threads(max_threads)
            {
                #pragma omp single
                while (itemset_x != itemsets.end())
                {

                    
                    #pragma omp task firstprivate(itemset_x), shared(temp)
                    {
                        auto itemset_y = itemset_x;
                        itemset_y++;
                        while (itemset_y != itemsets.end())
                        {
                            unsigned int *buf = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);
                            unsigned int pow_count = 0;

                            bool inserted = false, next = true;
                            __m128i *p_buf = (__m128i *)buf, *p_seti = (__m128i *)(*itemset_x), *p_setj = (__m128i *)(*itemset_y);
                            for (unsigned int mask = 0; next && pow_count <= 1 && mask < MASK_SIZE / 128; mask++)
                            {
                                unsigned int check_result[8];
                                __m128i m_seti = _mm_load_si128(p_seti);
                                __m128i or_result = _mm_or_si128(m_seti, _mm_load_si128(p_setj));
                                __m128i xor_result = _mm_xor_si128(m_seti, or_result);
                                _mm_storeu_si128((__m128i *)check_result, xor_result);
                                _mm_storeu_si128((__m128i *)(check_result + 4), _mm_and_si128(xor_result, _mm_sub_epi32(xor_result, _mm_set1_epi32(1))));
                                for (unsigned int b = 0; b < 4; b++)
                                {
                                    if (check_result[b])
                                        if (!(check_result[b + 4]))
                                            pow_count++;
                                        else
                                        {
                                            pow_count = 2;
                                            next = false;
                                        }
                                }
                                _mm_store_si128(p_buf, or_result);
                                ++p_buf;
                                ++p_seti;
                                ++p_setj;
                            }

                            if (pow_count == 1)
                            {
                                #pragma omp critical(merge_write)
                                    {inserted = temp.insert(buf).second;}
                            }
                            if (!inserted)
                                _mm_free(buf);
                            itemset_y++;
                            
                        }
                    }
                    ++itemset_x;
                }
                #pragma omp taskwait
                #pragma omp single
                {
                    itemsets.swap(temp);
                    for (auto set : temp)
                        #pragma omp task firstprivate(set)
                        {_mm_free(set);}
                }
                #pragma omp taskwait
            }
        }

        // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
    }
};


class SyncAprioriSSE : public SSE<VectorSSE>
{
    unsigned int cache_regulator;

    void singles_merge(const double support, const int max_threads)
    {
        cache_regulator = CACHE_SIZE / (MASK_SIZE/8); 
        #pragma omp parallel for schedule(dynamic) num_threads(max_threads)
        for (unsigned int i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned int j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        unsigned int *buf = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);
                        initialize_array(buf);
                        buf[i / 32] ^= static_cast<unsigned int>(std::exp2(31 - (i % 32) ));
                        buf[j / 32] ^= static_cast<unsigned int>(std::exp2(31 - (j % 32) ));
                        #pragma omp critical(single_write)
                        {
                            itemsets.push_back(buf);
                        }
                    }
                }
        }
    }

    void map(const unsigned int k, const int max_threads)
    {
        occurrencies = std::vector<unsigned int>(itemsets.size(), 0);
        if (transactions.size() < itemsets.size())
        {
            //for every itemset
            #pragma omp parallel num_threads(max_threads)
            for (unsigned int set = 0; set < itemsets.size(); set++)
            {
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned int tx = 0; tx < transactions.size(); tx++)
                {
                    unsigned int masks_number = masks[tx] / 128;
                    bool found = true;
                    __m128i *p_set = (__m128i *)itemsets[set], *p_tx = (__m128i *)transactions[tx];
                    for (unsigned int i = 0; i < MASK_SIZE / 128; i++)
                    {

                        __m128i set_128 = _mm_load_si128(p_set);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(p_tx)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                        ++p_set;
                        ++p_tx;
                    }
                    if (found)
                        #pragma omp atomic
                        occurrencies[set]++;
                }
                if(set % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
        else
        {
            //for every itemset
            #pragma omp parallel num_threads(max_threads)
            for (unsigned int tx = 0; tx < transactions.size(); tx++)
            {
                unsigned int masks_number = masks[tx] / 128;
                //for every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned int set = 0; set < itemsets.size(); set++)
                {
                    bool found = true;
                    __m128i *p_set = (__m128i *)itemsets[set], *p_tx = (__m128i *)transactions[tx];
                    for (unsigned int i = 0; i < MASK_SIZE / 128; i++)
                    {

                        __m128i set_128 = _mm_load_si128(p_set);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(p_tx)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                        ++p_set;
                        ++p_tx;
                    }
                    if (found)
                        #pragma omp atomic
                        occurrencies[set]++;
                }
                if (tx % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
    }

    void prune(const double support, const int max_threads )
    {   
        unsigned int size = transactions.size();
        unsigned int tail = itemsets.size() - 1;
        for(int i = tail; i>=0; i--)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) < support)
            {   
                _mm_free(itemsets[i]);
                itemsets[i] = itemsets[tail];
                --tail;
            }
        }
        itemsets.resize(tail+1);
        // #pragma omp parallel num_threads(max_threads)
        // #pragma single
        // #pragma omp task
        // {        
        //     store_itemsets("sse_output.dat");
        // }
    }

    void merge(const unsigned int k, const double support, const int max_threads)
    {
        // prune(support);
        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<unsigned int *, Compare> temp;
            std::vector<unsigned int *> v_temp;
            unsigned int size = transactions.size();
            unsigned int itemsets_size = itemsets.size();

            //for every itemset, try to unite it with another in the itemsets vector
            #pragma omp parallel num_threads(max_threads)
            {
                #pragma omp for schedule(dynamic)
                for (unsigned int i = 0; i < itemsets_size - 1; i++)
                {
                    if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                    {
                        for (unsigned int j = i + 1; j < itemsets_size; j++)
                        {

                            if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                            {
                                unsigned int *buf = (unsigned int *)_mm_malloc(MASK_SIZE / 8, 16);
                                unsigned int pow_count = 0;

                                bool inserted = false, next = true;
                                __m128i *p_buf = (__m128i *)buf, *p_seti = (__m128i *)itemsets[i], *p_setj = (__m128i *)itemsets[j];
                                for (unsigned int mask = 0; next && pow_count <= 1 && mask < MASK_SIZE / 128; mask++)
                                {
                                    unsigned int check_result[8];
                                    __m128i m_seti = _mm_load_si128(p_seti);
                                    __m128i or_result = _mm_or_si128(m_seti, _mm_load_si128(p_setj));
                                    __m128i xor_result = _mm_xor_si128(m_seti, or_result);
                                    _mm_storeu_si128((__m128i *)check_result, xor_result);
                                    _mm_storeu_si128((__m128i *)(check_result + 4), _mm_and_si128(xor_result, _mm_sub_epi32(xor_result, _mm_set1_epi32(1))));
                                    for (unsigned int b = 0; b < 4; b++)
                                    {
                                        if (check_result[b])
                                            if (!(check_result[b + 4]))
                                                pow_count++;
                                            else
                                            {
                                                pow_count = 2;
                                                next = false;
                                            }
                                    }
                                    _mm_store_si128(p_buf, or_result);
                                    ++p_buf;
                                    ++p_seti;
                                    ++p_setj;
                                }

                                if (pow_count == 1)
                                {
                                    #pragma omp critical(merge_write)
                                    if (inserted = temp.insert(buf).second)
                                        v_temp.push_back(buf);
                                }
                                if (!inserted)
                                    _mm_free(buf);
                            }
                        }
                    }
                }
                
                
                #pragma omp for schedule(static)
                for (unsigned int i = 0; i < itemsets.size(); i++)
                    _mm_free(itemsets[i]);
            }
            itemsets.swap(v_temp);
        }

        // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
    }

};