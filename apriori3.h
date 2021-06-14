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

std::vector<std::string> union_vec(
    const std::vector<std::string>& v1,
    const std::vector<std::string>& v2)
{
    // Declare a resultant union vector
    std::vector<std::string> union_vector;

    // Iterate through the input vector v1
    for (auto it = v1.begin(); it != v1.end(); it++)
        // For each element in v1, find an element
        // with the same value in the input vector v2
        if (std::find_if(v2.begin(), v2.end(), \
            [&](const std::string& str) { return str == *it; }) != v2.end())
                // If value of the current element in v1 is found in v2,
                // then add this value to the resultant union vector
                union_vector.push_back(*it);

    // Return the vector of union elements
    return union_vector;
}

std::vector<std::string> intersect_vec(
    const std::vector<std::string>& v1,
    const std::vector<std::string>& v2)
{
    // Declare a resultant intersection vector
    std::vector<std::string> intersect_vector;
    // Find a union of two vectors v1 and v2
    std::vector<std::string> union_vector = union_vec(v1, v2);

    // Declare vector v and fill it with element of vector v1
    std::vector<std::string> v(v1);
    // Insert elements of vector v2 into vector v
    v.insert(v.end(), v2.begin(), v2.end());

    // Sort elements of the vector v
    std::sort(v.begin(), v.end());
    // Prune all duplicate elements in vector v
    v.erase(std::unique(v.begin(), v.end()), v.end());

    // Iterate through the vector v
    for (auto it = v.begin(); it != v.end(); it++)
        // For each element in v, find an element
        // with the same value in the union vector
        if (std::find_if(union_vector.begin(), union_vector.end(),
            [&](const std::string& str) { return str == *it; }) == union_vector.end())
                // If value of the current element in v is found
                // in the union vector then add this value to
                // the resultant intersection vector
                intersect_vector.push_back(*it);

    // Return the intersection vector
    return intersect_vector;
}

std::vector<std::string> append_vec(
    const std::vector<std::string>& v1,
    const std::vector<std::string>& v2)
{
    // Declare vector v and fill it with element of vector v1
    std::vector<std::string> v(v1);
    // Insert elements of vector v2 into vector v
    v.insert(v.end(), v2.begin(), v2.end());
    // Sort elements of the vector v
    std::sort(v.begin(), v.end());
    // Prune all duplicate elements in vector v
    v.erase(std::unique(v.begin(), v.end()), v.end());

    // Return the resultant vector v
    return v;
}

bool compare_vec(
    const std::vector<std::string>& v1,
    const std::vector<std::string>& v2)
{
    // Declare vector v and assign it vector v2
    std::vector<std::string> v = v2;
    // Find union vector based on values from v1 and v2
    std::vector<std::string> v_u = union_vec(v1, v2);
    // If the union vector is not empty, return "true", otherwise "false"
    return (v_u.size() != 0) ? true : false;
}

bool contains(const std::vector<std::string>& v,
    const std::string& value)
{
    // Find the value in vector v
    // Return "true" if the value is found, otherwise "false"
    return std::find_if(v.begin(), v.end(),
        [&](const std::string& s) { return s == value; }) != v.end();
}

bool contains_vec(const std::vector<std::string>& v1,
    const std::vector<std::string>& v2)
{
    bool found = false;
    // Iterate through the vector v2,
    // until an element contained in both v1 and v2 is found
    for (auto it = v2.begin(); it != v2.end() && !found; it++)
        // For each element of v2, perform a check
        // if the current element of v2 is found in v1
        // As the result, assign found to "true" if the current
        // element is found and "false" unless otherwise
        found = contains(v1, *it) ? true : false;

    // Return the value of found variable
    return found;
}

std::vector<std::string> get_items(
    const std::vector<std::vector<std::string>>& T)
{
    // Declare a vector of items
    std::vector<std::string> items;

    // Iterate through the vector of transactions T
    for (auto it = T.begin(); it != T.end(); it++)
        // Insert each element from current
        // transaction itemset into the following items vector
        items.insert(items.end(), it->begin(), it->end());

    // Sort vector of items
    std::sort(items.begin(), items.end());
    // Prune all duplicate items
    items.erase(std::unique(items.begin(), items.end()), items.end());

    return items;
}


///////////////////////////////////////////////////////////////////////////
class Apriori {
    public:
    size_t lattice_cardinality=0;
    std::vector<std::vector<std::string>> transactions;
    std::vector<std::set<std::string>> itemsets;
    std::unordered_set<std::string> single_items;
    std::unordered_map<std::string, size_t> hashmap;
    



    inline void print_single_items (){
        for (auto i : single_items)
            std::cout <<"Element: " << i << " - Cardinality: "<< hashmap[i] << "\n";
        std::cout << "Total size: " << single_items.size() << "\n";
    }

    inline void print_items () {
        for (auto set : itemsets){
            for (auto str : set){
                std::cout << str <<"-";
            }
            std::cout << "\n";
        }
    }

    inline void update_single_structures(std::string token){
        single_items.insert(token);
        hashmap[token]++;       
    }

    void singles_prune(float support){
        auto item = std::begin(single_items);
        while (item != std::end(single_items)) {
            if (hashmap[*item] < support)
                item = single_items.erase(item);
            else
                ++item;
        }
        lattice_cardinality = single_items.size();
    }

    void singles_merge(){
        auto item_x = std::begin (single_items);
        auto item_y = item_x;
        while (item_x != std::end(single_items)){
            ++item_y;
            while (item_y != std::end(single_items)){
                itemsets.push_back({*item_x, *item_y});
                item_y++;
            }
            ++item_x;
            item_y = item_x;
        }
    }

    void map (float support){
        // std::cout << "ENTER MAP\n";
        hashmap.clear();
        //MAP
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
    }
    void prune(float support) {
        //PRUNE
        // std::cout << "ENTER PRUNE\n";
        std::cout << "ITEMSETS BEFORE PRUNING: " << itemsets.size() << "\n";
        auto itemset = std::begin(itemsets);
        while (itemset != std::end(itemsets)) {
            std::string buff;
            auto item = std::begin(*itemset);
            while (item != std::end(*itemset)){
                buff += (*item);
                item++;
            }
            if (hashmap[buff] < support)
                itemset = itemsets.erase(itemset);
            else
                ++itemset;
        }
        std::cout << "ITEMSETS AFTER PRUNING: " << itemsets.size() << "\n";
        // std::cout << "EXIT PRUNE\n";
    }

    void merge (int k){
        // std::cout << "ENTER MERGE\n";
        std::set<std::set<std::string>> temp;
        for (int itemset_x=0; itemset_x<itemsets.size()-1; itemset_x++){
            for (int itemset_y=itemset_x+1; itemset_y<itemsets.size(); itemset_y++){
                std::set<std::string> merged(itemsets[itemset_x]);
                merged.insert(itemsets[itemset_y].begin(), itemsets[itemset_y].end());
                if (merged.size() == k)
                    temp.insert(merged);
            }
        }
        itemsets.assign( temp.begin(), temp.end() );
        // std::cout << "EXIT MERGE\n";
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

            while(!getline(ifs, doc_buffer).eof()){
                std::vector<std::string> line_buffer;
                std::istringstream iss(doc_buffer);
                std::string token;
                while (std::getline(iss, token, ' ')){
                    line_buffer.push_back(token);
                    update_single_structures(token);                    
                }
                transactions.push_back(line_buffer);
            }
    }
    void run (const std::string input_file, float support){
        int k=3;
        read_data(input_file);
        singles_prune(support);
        singles_merge();
        print_single_items();
        while (!itemsets.empty()){     
            map(support);
            prune(support);
            merge(k);
            print_items();
            k++;
        }
    }
};






/* using APR_ENTITY = std::string;

std::vector<APR_ENTITY>
    read_sl_csv_file(const std::string filename)
{
    // Create an file input stream
    std::ifstream ifs(filename, \
        std::ifstream::in | std::ifstream::ate);

    // Declare a vector of entities
    std::vector<APR_ENTITY> v_ent;

    // Get the file size
    std::size_t file_size = (std::size_t)ifs.tellg();
    // Allocate a buffer to store the data from file
    char* buf = (char*)malloc(++file_size);
    // Read a single-line data from file
    ifs.seekg(ifs.beg); ifs.getline(buf, ++file_size);

    // Declare an entity object
    APR_ENTITY apr_ent;
    std::memset((void*)& apr_ent, 0x00, sizeof(APR_ENTITY));

    const char* delim = ",";
    // Find first token in buffer
    char* token = std::strtok(buf, delim);
    // Extract all tokens separated by ',' character
    while (token != nullptr)
    {
        std::uint32_t ci = 0L;
        // Allocate buffer to store an entity data
        char* m_value_buf = (char*)malloc(256 * sizeof(char));
        // Extract the data from the entity formatted string
        if (sscanf(token, "%[^:]:%d", m_value_buf, &ci) > 0)
        {
            // Assign the values extracted to the specific variables of entity object
            apr_ent.m_value = std::string(m_value_buf); apr_ent.m_ci = ci;
            // Add the entity object to the vector of entities
            // and extract the next token from the string data
            v_ent.push_back(apr_ent); token = std::strtok(nullptr, delim);
        }
    }

    // Deallocate buffer and close the input file stream
    free(buf); ifs.close();

    // Return the resultant vector of entities
    return v_ent;
}

std::vector<std::vector<std::string>> \
    read_data_from_csv_file(const std::string filename)
{
    // Declare a vector of transactions
    std::vector<std::vector<std::string>> trans_v;
    // Create an input file stream
    std::ifstream ifs(filename, std::ifstream::in);
    // Allocate buffer for each line of file
    char* buf = (char*)malloc(sizeof(char) * 4096);
    // Perform a check if this is not the end of file
    while (ifs.eof() == false && ifs.peek() >= 0)
    {
        // Declare a tokens buffer
        std::vector<std::string> tokens_v;
        // Read the current line from file
        const char* delim = ","; ifs.getline(buf, 4096);
        // Extract the first token from the current string
        char* token = std::strtok(buf, delim);
        // Extract all tokens separated by ',' character
        while (token != nullptr)
        {
            // Add the current token (i.e. item) to the vector of tokens
            tokens_v.push_back(std::string(token));
            // Extract the next token
            token = std::strtok(nullptr, delim);
        }

        // Add all items in the current vector to the vector of transactions
        trans_v.push_back(tokens_v);
    }

    // Deallocate buffer and close the input file stream
    free(buf); ifs.close();

    // Return the vector of transactions
    return trans_v;
}

void init(const std::vector<std::string>& items,
    std::vector<ASSOC_R_CAND>& R)
{
    // Declare a rule candidate object
    ASSOC_R_CAND asr_item;
    std::memset((void*)& asr_item, 0x00, sizeof(ASSOC_R_CAND));
    // Iterate through the vector of items
    for (auto it = items.begin(); it != items.end(); it++) {
        // For each item, assign the value of confidence equal to .0f
        asr_item.m_conf = .0f;
        // Create a single-element vector and assign it to the itemset variable
        asr_item.m_itemset = utils::vec1d(*it);
        // Add the current rule candidate object to the vector of rule candidates
        R.push_back(asr_item);
    }
}

    std::double_t get_conf_rule(
    const std::vector<std::vector<std::string>>& T,
    std::vector<std::string>& rule)
{
    // Declare a count variable and assign it to 0
    std::double_t count = 0L;
    // Iterate through the vector of transactions
    for (auto tran_it = T.begin(); tran_it != T.end(); tran_it++)
        // For each current transaction itemset, perform a check
        // if the union vector of the current itemset and the itemset of
        // the input rule has the same size as the input rule itemset vector
        if (utils::union_vec(*tran_it, rule).size() == rule.size())
            // If so, increment the count variable by 1
            count++;

    // Return the value of count variable
    return count;
}

void get_conf(
    const std::vector<std::vector<std::string>>& T,
    std::vector<ASSOC_R_CAND>& R)
{
    // Declare a count variable and assign it to 0
    std::double_t count = 0L;
    // Iterate through the vector of rule candidates
    for (auto asr_it = R.begin(); asr_it != R.end(); asr_it++)
        // Compute the value of count for the current rule candidate
        if ((count = get_conf_rule(T, asr_it->m_itemset)) > 0)
            // If the value of count is greater than 0, then
            // divide this value by the value of the overall transactions count
            // assign it to the current rule candidate m_conf member variable
            asr_it->m_conf = count / (std::double_t)T.size();
}



template<typename Predicate>
std::vector<ASSOC_R_CAND> find_rules(std::vector<ASSOC_R_CAND>& R, Predicate pred)
{
    // Declare a vector of rule candidates iterator
    // Declare a resultant vector of rule candidates
    // that satisfy a condition specified as lambda expression
    auto it = R.begin(); std::vector<ASSOC_R_CAND> ret_v;
    // Find each item in the vector of rule candidates that
    // satisfies a certain condition, and add it to the resultant
    // vector of rule candidates found
    while ((it = std::find_if(it, R.end(), pred)) != R.end())
        ret_v.push_back(*it++);

    // Return the resultant vector of rule candidates
    return ret_v;
}

bool contains_rule(const std::vector<ASSOC_R_CAND>& R,
    std::vector<std::string>& rule)
{
    bool found = false;
    // Declare a rule candidates vector iterator
    auto asr_it = R.begin();
    // Perform a search to find the first occurrence
    // of a rule candidate object contained in the
    // vector of rule candidates R
    while (asr_it != R.end() && !found)
        // For each rule candidate object, perform a check
        // if a union vector of both the current rule candidate itemset
        // and the input rule candidate itemset has the same size as the
        // input rule itemset. If so, assign the value of "true" to the found
        // variable, and assign "false", unless otherwise
        found = utils::union_vec((asr_it++)->m_itemset,
            rule).size() == rule.size() ? true : false;

    // Return the value of found variable
    return found;
}

std::vector<APR_RULE> \
		apriori(const std::vector<APR_ENTITY>& C,
			const std::vector<std::vector<std::string>>& T,	bool thorough = false)
	{
		// Declare a mutex lock object
		omp_lock_t s_lock;
		// Init the mutex lock object
		omp_init_lock(&s_lock);

		// Retrieve the vector of items
		std::vector<std::string> \
			items = utils::get_items(T);

		// Assign the value of the rooted class
		std::string rooted = C[0].m_value;

		// Declare a vector of rule candidates
		std::vector<ARL::ASSOC_R_CAND> rules_cand_v;
		// Init the vector of rule candidates and compute 
		// the confidence value for each specific rule candidate
		ARL::init(items, rules_cand_v); ARL::get_conf(T, rules_cand_v);

		// For each itemset of size k, generate rule candidates
		// and estimate the value of confidence of each new rule
		for (std::uint64_t k = 2; k <= items.size(); k++)
		{
			// Find all rule candidates which itemset size is equal to (k-1)
			std::vector<ARL::ASSOC_R_CAND> k_rule_cand_v = \
				ARL::find_rules(rules_cand_v, \
					[&](const ARL::ASSOC_R_CAND asr) { \
					return asr.m_itemset.size() == (k - 1); });

			// Declare a vector of new rule candidates
			std::vector<ARL::ASSOC_R_CAND> rules_cand_new_v;

			// Get the value of rule candidates vector size
			std::size_t k_rule_cand_v_size = k_rule_cand_v.size();

			// Generate pairs of new rule candidates (Note: This code is executed in parallel)
			#pragma omp parallel for collapse(2) schedule(guided, 4)
			for (std::uint64_t i = 0; i < k_rule_cand_v_size; i++)
				for (std::uint64_t j = 0; j < k_rule_cand_v_size; j++)
				{
					// Perform a check if the value of j is greater than the value of i
					// and we're not generating redunant equivalent rule candidates
					if (j <= i) continue;
					// For each pair of the rule candidates, perform a check
					// if the rooted item is not contain in either first or second
					// rule candidate itemset
					if ((utils::contains(k_rule_cand_v[i].m_itemset, rooted)) ||
						(utils::contains(k_rule_cand_v[j].m_itemset, rooted))) continue;

					// Perform a check if the value of confidence for both the first and second
					// rule candidate is greater than the specific confidence threshold value, or
					// if we're processing the 1 or 2-itemset rule candidates
					if ((k_rule_cand_v[i].m_conf > min_conf) &&
						(k_rule_cand_v[j].m_conf > min_conf) || ((k - 1) <= 2))
					{
						// Spawn a parallel task
						#pragma omp task
						{
							// If both rule candidates satisfies the condition above,
							// then perform a check if the first rule candidate does not
							// contain the second rule candidate. 
							if ((!utils::compare_vec(k_rule_cand_v[i].m_itemset,
								k_rule_cand_v[j].m_itemset)) || ((k - 1) <= 2))
							{
								// Declare a new rule candidate object
								ARL::ASSOC_R_CAND asr_rule;
								// If the condition above is true or we're still processing
								// the 1 or 2-itemset candidates, then merge both first and 
								// second rule candidates itemset
								asr_rule.m_itemset = utils::append_vec(\
									k_rule_cand_v[i].m_itemset, k_rule_cand_v[j].m_itemset);

								// Compute the value of support for the first candidate
								std::double_t rule_supp = \
									ARL::get_conf_rule(T, k_rule_cand_v[i].m_itemset);
								// Compute the value of support for the new rule candidate
								std::double_t rule_new_supp = \
									ARL::get_conf_rule(T, asr_rule.m_itemset);

								// Compute the value of confidence for the current new rule candidate
								asr_rule.m_conf = rule_new_supp / rule_supp;

								// Synchronize the following code by setting the mutex lock
								omp_set_lock(&s_lock);

								// Perform a check if the new rule candidate is not in the
								// vector of new rule candidates and the value of confidence
								// exceeds the specified confidence threshold
								if ((!ARL::contains_rule(rules_cand_new_v, \
									asr_rule.m_itemset)) && (asr_rule.m_conf > min_conf))
										// If so, add the rule candidate to the vector
										// of new rule candidates
										rules_cand_new_v.push_back(asr_rule);

								// Unset the mutex lock
								omp_unset_lock(&s_lock);
							}
						}
					}
				}

			// Insert the vector of new candidates produced into
			// the rule candidates vector
			rules_cand_v.insert(rules_cand_v.end(), \
				rules_cand_new_v.begin(), rules_cand_new_v.end());
		}

		// Declare an output vector
		std::vector<APR_RULE> output_v;
		
		// Iterate through the vector of classes 
		// (Note: This code is executed in parallel)
		#pragma omp parallel for schedule(dynamic, 4)
		for (auto c_it = C.begin() + 1; c_it != C.end(); c_it++)
		{
			// For each class value, find all rule candidates produced,
			// for which value of confidence exceeds the threshold specified
			// and an itemset each rule candidate contains a class value
			std::vector<ARL::ASSOC_R_CAND> rules_assoc_v = \
				ARL::find_rules(rules_cand_v, [&](const ARL::ASSOC_R_CAND& asr) { \
					return asr.m_conf > min_conf && \
						utils::contains(asr.m_itemset, c_it->m_value);
					});

			if (rules_assoc_v.size() > 0)
			{
				// Declare a vector of hold string values of
				// items from each rule candidate itemset that 
				// satisfies the condition above
				std::vector<std::string> \
					merged(rules_assoc_v.begin()->m_itemset);

				// Iterate through the resultant vector of rule 
				// candidates (Note: This code is executed in parallel)
				#pragma omp parallel for schedule(dynamic, 4)
				for (auto rl_it = rules_assoc_v.begin() + 1; rl_it != rules_assoc_v.end(); rl_it++)
					// Append all items in the itemset of each rule candidate
					// to the vector of items for the current class value
					merged = utils::append_vec(merged, rl_it->m_itemset);

				// Declare an association rule object
				APR_RULE apr_rule;
				std::memset((void*)& apr_rule, 0x00, sizeof(APR_RULE));
				// Assign the m_cs variable to a class name value
				apr_rule.m_cs = std::string(c_it->m_value);
				// Assign the m_itemset variable to the vector of merged items
				apr_rule.m_itemset = std::vector<std::string>(merged);

				// Remove an item, which value is equal to the value of class name
				apr_rule.m_itemset.erase(std::remove_if(\
					apr_rule.m_itemset.begin(), apr_rule.m_itemset.end(), \
						[&](const std::string s) { return s == c_it->m_value; }));

				// Add the rule to the output vector
				output_v.push_back(apr_rule);
			}
		}

		// Return the output vector
		return output_v;
	}
}

 */
