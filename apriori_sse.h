#include <vector>
#include <math.h>
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
#include <omp.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

using IndexType = unsigned int;

class SSE_Apriori
{

    std::vector<IndexType *> transactions;
    std::vector<IndexType *> itemsets;
    std::vector<std::string> single_items;
    std::vector<IndexType> occurrencies;
    IndexType mask_size = 128;

    inline void initialize_array(IndexType *buf)
    {
        for (IndexType m = 0; m < mask_size / 128; m++)
        {
            _mm_store_si128((__m128i *)(buf + (m * 4)), _mm_setzero_si128());
        }
    }


    inline void print_items()
    {
        for (IndexType is = 0; is < itemsets.size(); is++)
        {
            for (IndexType i = 0; i < mask_size / 32; i++)
            {
                std::bitset<32> x(itemsets[is][i]);
                std::cout << x << " ";
            }
            std::cout << "\n";
        }
    }

    inline void print_transactions()
    {
        for (IndexType is = 0; is < transactions.size(); is++)
        {
            for (IndexType i = 0; i < mask_size / 32; i++)
            {
                std::bitset<32> x(transactions[is][i]);
                std::cout << x << " ";
            }
            std::cout << "\n";
        }
    }

    void read_data(const std::string input_file, bool parallel=1)
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
            IndexType *buf = (IndexType *)_mm_malloc(mask_size / 8, 16);
            initialize_array(buf);
            for (IndexType m = 0; m < mask_size / 128; m++)
                transactions.push_back(buf);
            std::istringstream iss(doc_buffer);
            std::string token;
            while (std::getline(iss, token, ' '))
            {
                if (!isspace(token[0]))
                {
                    bool clear = 1;
                    for (IndexType i = 0; i < single_items.size(); i++)
                    {
                        if (single_items[i] == token)
                        {
                            clear = 0;
                            buf[i / 32] ^= static_cast<IndexType>(std::pow(2, (32 - (i % 32) - 1)));
                            if (occurrencies.size() < i)
                                occurrencies.resize(i + 1);
                            occurrencies[i]++;
                            // i=single_items.size();
                            break;
                        }
                    }
                    if (clear)
                    {
                        single_items.push_back(std::move(token));
                        occurrencies.push_back(1);
                        IndexType s_size = single_items.size();
                        if (single_items.size() > mask_size)
                        {
                            mask_size += 128;
                            IndexType *buf2 = (IndexType *)_mm_malloc(mask_size / 8, 16);
                            initialize_array(buf2);

                            for (IndexType i = 0; i < mask_size / 128; i++)
                            {
                                __m128i m = _mm_load_si128((__m128i *)buf + (i * 4));
                                _mm_store_si128((__m128i *)(buf2 + (i * 4)), m);
                            }
                            _mm_free(buf);
                            buf = buf2;
                        }
                        buf[(s_size - 1) / 32] ^= static_cast<IndexType>(std::pow(2, (32 - ((s_size - 1) % 32) - 1)));
                    }
                }
            }
        }
        if (mask_size != 128)
            #pragma omp parallel for if (parallel)
            for (IndexType i = 0; i < transactions.size(); i++)
            {
                IndexType *buf2 = (IndexType *)_mm_malloc(mask_size / 8, 16);
                initialize_array(buf2);
                for (IndexType j = 0; j < mask_size / 128; j++)
                {
                    __m128i m = _mm_load_si128((__m128i *)transactions[i] + (j * 4));
                    _mm_store_si128((__m128i *)(buf2 + (j * 4)), m);
                }
                _mm_free(transactions[i]);
                transactions[i] = buf2;
            }
    }

    void singles_merge(double support, bool parallel=1)
    {
        #pragma omp parallel for if (parallel) schedule(dynamic)
        for (IndexType i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (IndexType j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        IndexType *buf = (IndexType *)_mm_malloc(mask_size / 8, 16);
                        initialize_array(buf);
                        buf[i / 32] ^= static_cast<IndexType>(std::pow(2, (32 - (i % 32) - 1)));
                        buf[j / 32] ^= static_cast<IndexType>(std::pow(2, (32 - (j % 32) - 1)));
                        #pragma omp critical (single_write)
                        {itemsets.push_back(buf);}
                    }
                }
        }
    }

    void map(IndexType k, bool parallel=1)
    {

        occurrencies.resize(itemsets.size());
        //for every itemset
        #pragma omp parallel for if(parallel)
        for (IndexType set = 0; set < itemsets.size(); set++)
        {
            occurrencies[set] = 0;
            //for every transaction
            for (IndexType tx = 0; tx < transactions.size(); tx++)
            {
                IndexType found = 0;
                for (IndexType i = 0; i < mask_size / 128; i++)
                {
                    __m128i set_128 = _mm_load_si128((__m128i *)itemsets[set] + (i * 4));
                    __m128i xor_result = _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128((__m128i *)transactions[tx] + (i * 4))), set_128);
                    found += _mm_testz_si128(xor_result, xor_result);
                    if (found == i)
                        break;
                }
                if (found == mask_size / 128)
                    occurrencies[set]++;
            }
        }

    }



    inline void countSetBits(IndexType *buf, IndexType &bitcount)
    {

        __m128i sse_bitcount = _mm_set1_epi32(0); // accumulate 4 partial counters
        const __m128i mask = _mm_set1_epi32(0x00000001);
        __m128i *p_sse = (__m128i *)buf;
        for (IndexType i = 0; i < mask_size / 32; i += 4)
        { // 4 bytes at once
            __m128i copy = _mm_load_si128(p_sse);
            for (IndexType j = 0; j < 32; j++)
            {
                __m128i result = _mm_and_si128(copy, mask);
                sse_bitcount = _mm_add_epi32(sse_bitcount, result);
                copy = _mm_srli_epi32(copy, 1);
            }
            p_sse++;
        }

        // sum the four partial results
        unsigned int partial[4];
        _mm_storeu_si128((__m128i *)partial, sse_bitcount);
        bitcount += partial[0] + partial[1] + partial[2] + partial[3];
    }

    void merge(unsigned int k, double support, bool parallel=1)
    {

        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            std::set<std::vector<IndexType>> temp;
            std::vector<IndexType *> v_temp;
            IndexType size = transactions.size();
            IndexType itemsets_size = itemsets.size();
            
            //for every itemset, try to unite it with another in the itemsets vector
            #pragma omp parallel for if(parallel) schedule(dynamic)
            for (IndexType i = 0; i < itemsets_size - 1; i++)
            {
                // std::cout << static_cast<double>(occurrencies[i]) << " ";
                if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                {
                    // std::cout << i << " ";
                    for (IndexType j = i + 1; j < itemsets_size; j++)
                    {

                        if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                        {
                            IndexType *buf = (IndexType *)_mm_malloc(mask_size / 8, 16);
                            initialize_array(buf);
                            IndexType bit_count = 0, aggregator = 0;

                            for (IndexType mask = 0; mask < mask_size / 32; mask += 4)
                            {
                                _mm_store_si128((__m128i *)(buf + mask), _mm_or_si128(_mm_load_si128((__m128i *)itemsets[i] + mask), _mm_load_si128((__m128i *)itemsets[j] + mask)));
                                countSetBits(buf + mask, bit_count);
                            }
                            if (bit_count == k)
                            {
                                
                                #pragma omp critical (merge_write)
                                if (temp.insert(std::vector<IndexType>(buf, buf + mask_size / 32)).second)
                                
                                    v_temp.push_back(buf);
                                
                            }
                            else
                                _mm_free(buf);
                        }
                    }
                }
            }
            itemsets.swap(v_temp);
            #pragma omp parallel for if (parallel)
            for (IndexType i=0; i<v_temp.size(); i++)
                _mm_free(v_temp[i]);
        }

        // std::cout << "ITEMSETS SIZE: " << itemsets.size() <<  "\n";
    }

public:
    void run(const std::string input_file, double support, bool parallel=1)
    {
        IndexType k = 2;
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
