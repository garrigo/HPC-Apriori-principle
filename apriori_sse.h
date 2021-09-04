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

//Number of masks of 128 bits used when creating or comparing itemsets
static unsigned MASK_SIZE;

//Comparator of 2 itemsets using 2 SSE cmpgt for every mask
struct Compare
{
    bool operator()(unsigned * __restrict__ a, unsigned * __restrict__ b) const
    {
        /*Array of results of comparisons
        Example: if compare_result[0] > 0, first 1/4 of the mask of a is lesser than first 1/4 of mask of b
        If compare_result[4] > 0 we have the opposite situation
        If neither are true, the 1/4 of masks are equal*/
        unsigned compare_result[8];

        //For every mask of the itemset, starting from the first one
        for (unsigned mask = 0; mask < MASK_SIZE; mask++)
        {
            //Load a and b masks number [mask]
            __m128i m_a = _mm_load_si128(((__m128i *)a)+mask);
            __m128i m_b = _mm_load_si128(((__m128i *)b)+mask);

            /*Store the comparisons of the two loop's masks in the two halves of the array
            The second one has the order of arguments inverted (first checks b>a, second a>b)*/
            _mm_storeu_si128((__m128i *)compare_result, _mm_cmpgt_epi32(m_b, m_a));
            _mm_storeu_si128((__m128i *)(compare_result + 4), _mm_cmpgt_epi32(m_a, m_b));
            for (unsigned b = 0; b < 4; b++)
            {
                //If the first mask's leftmost values are smaller than the second ones, the first itemset is smaller
                if (compare_result[b])
                    return true;
                if (compare_result[b + 4])
                    return false;
            }
        }
        //The itemsets have all masks equal
        return false;
    }
};

//Data type of the itemset collection for the Sync version
struct VectorSSE {
    typedef std::vector<unsigned *> data_type;
};

//Data type of the itemset collection for the normal version
struct SetSSE {
    typedef std::set<unsigned *, Compare> data_type;
};

//Abstract class that includes all members that are shared by all SSE versions
template<class ItemsetType>
class SSE
{
protected:
    //Vector of unique items taken by the file
    std::vector<std::string> single_items;
    //Vector of masks of bits representing the single items present in a transaction
    std::vector<unsigned *> transactions;
    //Vector of mask sizes for each transaction
    std::vector<unsigned> masks;
    //Collection of itemsets generated in the k-step
    typename ItemsetType::data_type  itemsets;
    //Vector representing the occurrencies of every itemset in the whole transaction list
    std::vector<unsigned> occurrencies;

    //Methods to be implemented by every child class

    virtual void singles_merge(const double, const int) = 0;
    virtual void map(const double, const unsigned, const int) = 0;
    virtual void merge(const unsigned, const double, const int) = 0;
  
    //Output single items vector to file
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

    //Output current itemsets to file
    void store_itemsets(std::string filename)
    {
        std::ofstream ofs;
        ofs.open(filename, std::ios_base::app);
        if (ofs.is_open())
        {        
            for (auto& set : itemsets)
            {
                for (unsigned i = 0; i < MASK_SIZE*4; i++)
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

    //Set all bits of an array to 0 using SSE
    inline void initialize_array(unsigned *buf)
    {
        for (unsigned m = 0; m < MASK_SIZE; m++)
            _mm_store_si128((__m128i *)(buf + (m * 4)), _mm_setzero_si128());
    }

    //Read the input file and populate single_items, transactions and occurrencies
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
        #pragma omp parallel num_threads(max_threads)
        #pragma omp single
        while (!getline(ifs, doc_buffer).eof())
        //Create a task for every line of the input file
        #pragma omp task firstprivate(doc_buffer)
        {
            
            {
                std::istringstream iss(doc_buffer);
                std::string token;
                //Array of indices of tokens found w.r.t. single_items vector
                std::vector<unsigned> indices;
                //Max index value in the indices vector
                unsigned max = 0;
                //For every token in a line
                while (std::getline(iss, token, ' '))
                {
                    if (!isspace(token[0]))
                    {
                        bool clear = true;
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
                        //Save the position of the token (in the single_items vector) in the temporary vector of indices
                        indices.push_back(i);
                        if ( i > max)
                            max = i;
                    }
                }
                //Number of 128 bit masks to be allocated for the transaction
                unsigned mask = max/128+1;
                unsigned * buf = (unsigned * )_mm_malloc(mask*16, 16);
                //Set all bits to 0
                for (unsigned m = 0; m <= max / 128; m++)
                    _mm_store_si128((__m128i *)(buf + (m * 4)), _mm_setzero_si128());
                //Set to 1 every bit in the [index] position of the buf array       
                for (auto& i: indices)
                    buf[i / 32] ^= static_cast<unsigned>(std::exp2(31 - (i % 32) ));
                //Push the transaction and its mask size and modify the global MASK_SIZE if needed
                #pragma omp critical (mask_write)
                {
                    transactions.push_back(buf);
                    masks.push_back(mask);
                    if (MASK_SIZE < max)
                        MASK_SIZE = mask;
                }
            }
        }
        ifs.close();
    } 

public:
    void run(const std::string input_file, const double support, const int max_threads)
    {   
        
        unsigned k = 2;
        MASK_SIZE = 1;
        read_data(input_file, max_threads);
        singles_merge(support, max_threads);
        //While itemsets within the given threshold are being found
        while (!itemsets.empty())
        {
            map(support, k, max_threads);
            ++k;
            merge(k, support, max_threads);
        }
        //De-allocate from memory the arrays of transactions
        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for (unsigned i = 0; i < transactions.size(); i++)
            _mm_free(transactions[i]);
    }
};



class AprioriSSE : public SSE<SetSSE>
{
    //Create itemsets by merging couples of single items (for the k=2 pass)
    void singles_merge(const double support, const int max_threads)
    {
        #pragma omp parallel for num_threads(max_threads) schedule(dynamic, 1)
        for (unsigned i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        unsigned *buf = (unsigned *)_mm_malloc(MASK_SIZE * 16, 16);
                        initialize_array(buf);
                        buf[i / 32] ^= static_cast<unsigned>(std::exp2(31 - (i % 32) ));
                        buf[j / 32] ^= static_cast<unsigned>(std::exp2(31 - (j % 32) ));
                        #pragma omp critical(single_write)
                        {
                            itemsets.insert(buf);
                        }
                    }
                }
        }
      
        // store_single_items("sse_output.dat");
        // store_itemsets("sse_output.dat");
    }

    //For every itemset, search it in every transaction and if it's present increment the corresponding occurrencies' entry
    void map(const double support, const unsigned k, const int max_threads)
    {
        unsigned occ = 0;
        occurrencies.resize(itemsets.size());
        //For every itemset create a task
        #pragma omp parallel num_threads(max_threads)
        #pragma omp single
        for (const auto set : itemsets)
        {
            occurrencies[occ] = 0;
            //For every transaction
            #pragma omp task firstprivate(set, occ)
            {
                for (unsigned tx = 0; tx < transactions.size(); tx++)
                {
                    bool found = true;
                    unsigned masks_number = masks[tx];
                    __m128i* tx_128 = (__m128i*)transactions[tx];
                    /*Check if the 1 bits in the [i] mask of the itemset are all included in the [i] mask of the transaction
                    If the number of masks of the transaction is smaller than the itemset's one, test if the itemset's mask is 0*/
                    for (unsigned i = 0; i < MASK_SIZE; i++)
                    {
                        __m128i set_128 = _mm_load_si128((__m128i *)(set)+i);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(tx_128+i)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                    }
                    if (found)
                        occurrencies[occ]++;
                }
            }
            ++occ;
        }
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
                _mm_free(*item_set);
                item_set = itemsets.erase(item_set);
            }
            else
                ++item_set;
            ++occ;
        }
      
        // store_itemsets("sse_output.dat");

    }

    //Merge all k-1 pass itemsets to create a new set of itemsets having k single items
    void merge(const unsigned k, const double support, const int max_threads)
    {
        prune(support, max_threads);
        if (!itemsets.empty())
        {
            std::set<unsigned *, Compare> temp;
            auto itemset_x = itemsets.begin();
            
            #pragma omp parallel num_threads(max_threads)
            {
                //A single thread creates a task for every itemset
                #pragma omp single
                while (itemset_x != itemsets.end())
                {
                    //The task tries to unite the given itemset with all the itemsets following it
                    #pragma omp task firstprivate(itemset_x), shared(temp)
                    {
                        auto itemset_y = itemset_x;
                        itemset_y++;
                        while (itemset_y != itemsets.end())
                        {
                            //New array to contain the new itemset
                            unsigned *buf = (unsigned *)_mm_malloc(MASK_SIZE * 16, 16);
                            unsigned pow_count = 0;
                            bool inserted = false;
                            /*For every mask of the two itemsets compute the or (union) and to check if the new itemset has exactly k 1 bits:
                            1) Compute the xor between one original itemset and the or result and store it in the first half of the check_result array
                            2) Subtract 1 to every chunk of 32 bits of the xor result.
                            3) Compute the and operation between the xor result and the sub result and store it in the second half of the check_result array
                            This computation checks if the 128 bit mask resulting from the xor gives a number that is a power of 2 (only 1 bit set),
                            which means that in order for the new itemset to have k 1 bits, there must be just one bit set in all masks of the xor result. */
                            for (unsigned mask = 0; pow_count <= 1 && mask < MASK_SIZE; mask++)
                            {
                                unsigned check_result[8];
                                __m128i m_seti = _mm_load_si128(((__m128i*)(*itemset_x))+mask);
                                __m128i or_result = _mm_or_si128(m_seti, _mm_load_si128(((__m128i*)(*itemset_y))+mask));
                                __m128i xor_result = _mm_xor_si128(m_seti, or_result);
                                _mm_storeu_si128((__m128i *)check_result, xor_result);
                                _mm_storeu_si128((__m128i *)(check_result + 4), _mm_and_si128(xor_result, _mm_sub_epi32(xor_result, _mm_set1_epi32(1))));
                                for (unsigned b = 0; b < 4; b++)
                                {
                                    if (check_result[b])
                                        if (!(check_result[b + 4]))
                                            pow_count++;
                                        else
                                        {
                                            pow_count = 2;
                                            break;
                                        }
                                }
                                _mm_store_si128(((__m128i*)buf)+mask, or_result);
                            }
                            //If the xor masks produce just 1 bit, the itemset have k 1 bits and can be inserted
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
                //Create a task for every old itemset to free them
                #pragma omp single
                {
                    itemsets.swap(temp);
                    for (auto set : temp)
                        #pragma omp task firstprivate(set)
                        {_mm_free(set);}
                }
            }
        }

        // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
    }
};


class SyncAprioriSSE : public SSE<VectorSSE>
{
    //Custom value for the computation of cache regulator (cache aware model)
    static constexpr size_t CACHE_SIZE = 1048576;
    unsigned cache_regulator;

    //Create itemsets by merging couples of single items (for the k=2 pass)
    void singles_merge(const double support, const int max_threads)
    {
        //Number of lines to synchronize
        cache_regulator = CACHE_SIZE / (MASK_SIZE * 16); 
        
        #pragma omp parallel num_threads(max_threads)
        #pragma omp for schedule(dynamic, 1) 
        for (unsigned i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                // #pragma omp for schedule(static) nowait
                for (unsigned j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        unsigned *buf = (unsigned *)_mm_malloc(MASK_SIZE * 16, 16);
                        initialize_array(buf);
                        buf[i / 32] ^= static_cast<unsigned>(std::exp2(31 - (i % 32) ));
                        buf[j / 32] ^= static_cast<unsigned>(std::exp2(31 - (j % 32) ));
                        #pragma omp critical(single_write)
                        {
                            itemsets.push_back(buf);
                        }
                    }
                }
            //Every [cache_regulator] items, synchronize threads making the fastest ones wait 
            // if (i % cache_regulator == 0)
            // {
            //     #pragma omp barrier
            // }
        }
        // store_single_items("sse_output.dat");
        // store_itemsets("sse_output.dat");
    }

    //For every itemset, search it in every transaction and if it's present increment the corresponding occurrencies' entry
    void map(const double support, const unsigned k, const int max_threads)
    {
        occurrencies = std::vector<unsigned>(itemsets.size(), 0);
        if (transactions.size() < itemsets.size())
        {
            //For every itemset
            #pragma omp parallel num_threads(max_threads)
            for (unsigned set = 0; set < itemsets.size(); set++)
            {
                //For every transaction
                #pragma omp for schedule(static) nowait
                for (unsigned tx = 0; tx < transactions.size(); tx++)
                {
                    unsigned masks_number = masks[tx];
                    bool found = true;
                    /*Check if the 1 bits in the [i] mask of the itemset are all included in the [i] mask of the transaction
                    If the number of masks of the transaction is smaller than the itemset's one, test if the itemset's mask is 0*/                    
                    for (unsigned i = 0; i < MASK_SIZE; i++)
                    {
                        __m128i set_128 = _mm_load_si128(((__m128i *)itemsets[set])+i);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(((__m128i *)transactions[tx])+i)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                    }
                    if (found)
                        #pragma omp atomic
                        occurrencies[set]++;
                }
                //Synchronize threads every[cache_regulator] line
                if(set % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
        else
        {
            //For every transaction
            #pragma omp parallel num_threads(max_threads)
            for (unsigned tx = 0; tx < transactions.size(); tx++)
            {
                unsigned masks_number = masks[tx];
                //For every itemset
                #pragma omp for schedule(static) nowait
                for (unsigned set = 0; set < itemsets.size(); set++)
                {
                    bool found = true;
                    for (unsigned i = 0; i < MASK_SIZE; i++)
                    {
                    /*Check if the 1 bits in the [i] mask of the itemset are all included in the [i] mask of the transaction
                    If the number of masks of the transaction is smaller than the itemset's one, test if the itemset's mask is 0*/                    
                        __m128i set_128 = _mm_load_si128(((__m128i *)itemsets[set])+i);
                        __m128i xor_result = (masks_number > i) ? _mm_xor_si128(_mm_and_si128(set_128, _mm_load_si128(((__m128i *)transactions[tx])+i)), set_128) : set_128;
                        if (!(found = _mm_testz_si128(xor_result, xor_result)))
                            break;
                    }
                    if (found)
                        #pragma omp atomic
                        occurrencies[set]++;
                }
                //Synchronize threads every[cache_regulator] line
                if (tx % cache_regulator == 0)
                {
                    #pragma omp barrier
                }
            }
        }
    }

    //Prune all itemsets that do not reach the threshold
    void prune(const double support, const int max_threads )
    {   
        unsigned size = transactions.size();
        //Populate a new vector with the itemsets
        std::vector<unsigned *> temp;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(max_threads)
        for(int i = 0; i<itemsets.size(); i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
            {   
                #pragma omp critical (prune)
                {temp.push_back(itemsets[i]);}
            }
            else
                _mm_free(itemsets[i]);
        }
        itemsets.swap(temp);
  
        // store_itemsets("sse_output.dat");
    }

    void merge(const unsigned k, const double support, const int max_threads)
    {
        prune(support, max_threads);
        if (!itemsets.empty())
        {
            std::set<unsigned *, Compare> temp;
            std::vector<unsigned *> v_temp;
            unsigned size = transactions.size();
            unsigned itemsets_size = itemsets.size();

            //Try to unite every itemset with the following ones in the vector
            #pragma omp parallel num_threads(max_threads)
            {
                #pragma omp for schedule(dynamic, 1)
                for (unsigned i = 0; i < itemsets_size - 1; i++)
                {
                    // if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(size))) >= support)
                    {
                        // #pragma omp for schedule(static) nowait
                        for (unsigned j = i + 1; j < itemsets_size; j++)
                        {

                            // if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(size))) >= support)
                            {
                                //New array to contain the new itemset
                                unsigned *buf = (unsigned *)_mm_malloc(MASK_SIZE * 16, 16);
                                unsigned pow_count = 0;
                                bool inserted = false;
                                /*For every mask of the two itemsets compute the or (union) and to check if the new itemset has exactly k 1 bits:
                                1) Compute the xor between one original itemset and the or result and store it in the first half of the check_result array
                                2) Subtract 1 to every chunk of 32 bits of the xor result.
                                3) Compute the and operation between the xor result and the sub result and store it in the second half of the check_result array
                                This computation checks if the 128 bit mask resulting from the xor gives a number that is a power of 2 (only 1 bit set),
                                which means that in order for the new itemset to have k 1 bits, there must be just one bit set in all masks of the xor result. */
                                for (unsigned mask = 0; pow_count <= 1 && mask < MASK_SIZE; mask++)
                                {
                                    unsigned check_result[8];
                                    __m128i m_seti = _mm_load_si128((__m128i*)(itemsets[i])+mask);
                                    __m128i or_result = _mm_or_si128(m_seti, _mm_load_si128((__m128i*)(itemsets[j])+mask));
                                    __m128i xor_result = _mm_xor_si128(m_seti, or_result);
                                    _mm_storeu_si128((__m128i *)check_result, xor_result);
                                    _mm_storeu_si128((__m128i *)(check_result + 4), _mm_and_si128(xor_result, _mm_sub_epi32(xor_result, _mm_set1_epi32(1))));
                                    for (unsigned b = 0; b < 4; b++)
                                    {
                                        if (check_result[b])
                                            if (!(check_result[b + 4]))
                                                pow_count++;
                                            else
                                            {
                                                pow_count = 2;
                                                break;
                                            }
                                    }
                                    _mm_store_si128(((__m128i*)buf)+mask, or_result);
                                }
                                //If the xor masks produce just 1 bit, the itemset have k 1 bits and can be inserted
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
                    // if (i % cache_regulator == 0)
                    // {
                    //     #pragma omp barrier
                    // }
                }
                #pragma omp barrier
                //De-allocate the old itemsets
                #pragma omp for schedule(static)
                for (unsigned i = 0; i < itemsets.size(); i++)
                    _mm_free(itemsets[i]);
            }
            itemsets.swap(v_temp);
        }

        // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
    }

};