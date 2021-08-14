#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_set>
#include <omp.h>


class Apriori
{

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
    using U_VectorSet = std::unordered_set<std::vector<unsigned int>, VectorHash>;

    std::vector<std::vector<unsigned int>> transactions;
    U_VectorSet itemsets;
    std::vector<std::string> single_items;
    std::vector<unsigned int> occurrencies;

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
        unsigned int current_size = 0;
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
                            if (occurrencies.size() < i)
                                occurrencies.resize(i + 1);
                            occurrencies[i]++;
                            i = single_items.size();
                            break;
                        }
                    }
                    if (clear)
                    {
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

    void singles_prune(double support)
    {
    }

    void singles_merge(double support)
    {

        for (unsigned int i = 0; i < single_items.size() - 1; i++)
        {
            if ((static_cast<double>(occurrencies[i]) / (static_cast<double>(transactions.size()))) >= support)
                for (unsigned int j = i + 1; j < single_items.size(); j++)
                {
                    if ((static_cast<double>(occurrencies[j]) / (static_cast<double>(transactions.size()))) >= support)
                    {
                        itemsets.insert({i, j});
                    }
                }
        }
    }

    void map(unsigned int k)
    {
        occurrencies.resize(itemsets.size());
        //for every itemset
        unsigned int occ = 0;
        for (const auto set : itemsets)
        {
            occurrencies[occ] = 0;
            //for every transaction
            for (auto &tx : transactions)
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
                            // if ((*item)<(tx[tx_cursor])){
                            //     tx_cursor=tx.size();
                            // }
                            // else
                            if ((*item) == (tx[tx_cursor]))
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
            ++occ;
        }
    }

    void prune(double support)
    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
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
    }

    void merge(int k, double support)
    {

        if (!itemsets.empty())
        {
            // provisional set of set of string to modify the current itemsets vector with k+1 cardinality
            U_VectorSet temp;
            //for every itemset, try to unite it with another in the itemsets vector
            auto itemset_x = itemsets.begin();
            unsigned int size = transactions.size();
            while (itemset_x != itemsets.end())
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
                        temp.insert(merged);
                    itemset_y++;
                }

                ++itemset_x;
            }
            itemsets.swap(temp);
            // std::cout << "ITEMSETS SIZE: " << itemsets.size() << "\n";
        }
    }

public:
    void run(const std::string input_file, double support)
    {
        int k = 2;
        read_data(input_file);
        singles_merge(support);
        while (!itemsets.empty())
        {
            map(k);
            prune(support);
            ++k;
            merge(k, support);
        }
    }
};
