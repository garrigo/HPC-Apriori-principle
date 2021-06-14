#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <tuple>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using Transaction = std::vector<int>;

char DELIMITER = ' ';

inline long double round(long double value, int pos){
    long double temp;
    temp = value * pow( 10, pos );
    temp = floor( temp + 0.5 );
    temp *= pow( 10, -pos );
    return temp;
}

class Apriori {
private:
    int nowStep;
    long double minSupport;
    std::vector <Transaction> transactions;
    std::vector <Transaction > C, L;
    std::vector <std::vector <Transaction > > frequentSet;
    std::vector <tuple <Transaction, Transaction, long double, long double> > associationRules;
public:
    Apriori (std::vector<Transaction > _transactions, long double _minSupport) {
        nowStep=0;
        minSupport = _minSupport;
        for(auto&row:_transactions){
            sort(row.begin(), row.end());
            transactions.push_back(row);
        }
        frequentSet.push_back({{}});
        std::cout << "APRIORI INSTANTIATION ENDED\n";
    }
    
    std::vector <tuple <Transaction, Transaction, long double, long double> > getAssociationRules(){
        return associationRules;
    }
    
    void process() {
        while(true) {
            C = generateNextC();
            if(C.size()==0) break;
            nowStep++;
            
            L = generateL();
            frequentSet.push_back(L);
        }
        
        for(auto&stepItemSet:frequentSet) {
            for(auto&items:stepItemSet) {
                generateAssociationRule(items, {}, {}, 0);
            }
        }
    }
    
    void generateAssociationRule(Transaction items, Transaction X, Transaction Y, int index) {
        if(index == items.size()) {
            if(X.size()==0 || Y.size() == 0) return;
            long double XYsupport = getSupport(items);
            long double Xsupport = getSupport(X);
            
            if(Xsupport == 0) return;
            
            long double support = (long double)XYsupport;
            long double confidence = (long double)XYsupport/Xsupport*100.0;
            associationRules.push_back({X, Y, support, confidence});
            return;
        }
        
        X.push_back(items[index]);
        generateAssociationRule(items, X, Y, index+1);
        X.pop_back();
        Y.push_back(items[index]);
        generateAssociationRule(items, X, Y, index+1);
    }
    
    Transaction getElement(std::vector<Transaction > itemset) {
        Transaction element;
        set<int> s;
        for(auto&row:itemset) for(auto&col:row) s.insert(col);
        for(auto iter=s.begin(); iter != s.end(); iter++) element.push_back(*iter);
        return element;
    }
    
    std::vector<Transaction > generateNextC() {
        if(nowStep==0) {
            std::vector<Transaction > ret;
            Transaction element = getElement(transactions);
            for(auto&i:element) ret.push_back(Transaction(1, i));
            
            return ret;
        } else {
            return pruning(joining());
        }
    }
    
    std::vector<Transaction > joining () {
        std::vector<Transaction > ret;
        for(int i=0;i<L.size();i++){
            for(int j=i+1;j<L.size(); j++) {
                int k;
                for(k=0;k<nowStep-1; k++) {
                    if(L[i][k] != L[j][k]) break;
                }
                if(k == nowStep-1) {
                    Transaction tmp;
                    for(int k=0;k<nowStep-1; k++) {
                        tmp.push_back(L[i][k]);
                    }
                    int a = L[i][nowStep-1];
                    int b = L[j][nowStep-1];
                    if(a>b) swap(a,b);
                    tmp.push_back(a), tmp.push_back(b);
                    ret.push_back(tmp);
                }
            }
        }
        return ret;
    }
    
    std::vector<Transaction > pruning (std::vector<Transaction > joined) {
        std::vector<Transaction > ret;
        
        set<Transaction > lSet;
        for(auto&row:L) lSet.insert(row);
        
        for(auto&row:joined){
            int i;
            for(i=0;i<row.size();i++){
                Transaction tmp = row;
                tmp.erase(tmp.begin()+i);
                if(lSet.find(tmp) == lSet.end()) {
                    break;
                }
            }
            if(i==row.size()){
                ret.push_back(row);
            }
        }
        return ret;
    }
    
    long double getSupport(Transaction item) {
        int ret = 0;
        for(auto&row:transactions){
            int i, j;
            if(row.size() < item.size()) continue;
            for(i=0, j=0; i < row.size();i++) {
                if(j==item.size()) break;
                if(row[i] == item[j]) j++;
            }
            if(j==item.size()){
                ret++;
            }
        }
        return (long double)ret/transactions.size()*100.0;
    }
    
    std::vector<Transaction > generateL() {
        std::vector<Transaction > ret;
        for(auto&row:C){
            long double support = getSupport(row);
            if(round(support, 2) < minSupport) continue;
            ret.push_back(row);
        }
        return ret;
    }
};

class InputReader {
    private:
        ifstream input_file;
        std::vector<Transaction > transactions;
    public:
        InputReader(string filename) {
            input_file.open(filename);
            if(!input_file) {
                cout << "Input file could not be opened\n";
                exit(0);
            }
            string str_buffer;
            const char delimiter = ' ';
            while(!getline(input_file, str_buffer).eof()){
                Transaction int_buffer;
                int pre = 0;
                for(int i=0;i<str_buffer.size()-2;i++){
                    if(str_buffer[i] == DELIMITER) {
                        int num = atoi(str_buffer.substr(pre, i).c_str());
                        int_buffer.push_back(num);
                        pre = i+1;
                    }
                }
                int num = atoi(str_buffer.substr(pre, str_buffer.size()).c_str());
                int_buffer.push_back(num);
                
                // Extract the first token from the current string
                // int num;
                // char* token = std::strtok(&str_buffer[0], &delimiter);
                // int_buffer.push_back(atoi(token));       
                // while (token != nullptr){
                //     // Extract the next token
                //     token = std::strtok(nullptr, &delimiter);
                //     // Add the current token (i.e. item) to the vector of tokens
                //     int_buffer.push_back(atoi(token));
            
            
                // } 
                transactions.push_back(int_buffer);
            }
        }

        std::vector<Transaction > getTransactions() {
            // int r = 1;
            // for (auto& row : transactions){
                
            //     std::cout << r << "° Line: ";
            //     for (auto& item: row)
            //         std::cout  << item << " ";
            //     std::cout << "\n";
            //     r++;
            // }
            return transactions;
        }
};

class OutputPrinter {
private:
    ofstream fout;
    std::vector<tuple<Transaction, Transaction, long double, long double> > associationRules;
public:
    OutputPrinter(string filename, std::vector<tuple<Transaction, Transaction, long double, long double> > _associationRules) {
        fout.open(filename);
        if(!fout) {
            cout << "Ouput file could not be opened\n";
            exit(0);
        }
        associationRules = _associationRules;
        buildOutput();
    }
    
    void buildOutput() {
        for(auto&i:associationRules) {
            fout << vectorToString(get<0>(i)) << '\t';
            fout << vectorToString(get<1>(i)) << '\t';
            
            fout << fixed;
            fout.precision(2);
            fout << get<2>(i) << '\t';
            
            fout << fixed;
            fout.precision(2);
            fout << get<3>(i);
            
            fout << endl;
        }
    }
    
    string vectorToString(Transaction arr) {
        string ret = "{";
        for(int i=0;i<arr.size();i++){
            ret += to_string(arr[i]);
            if(i != arr.size()-1){
                ret += ",";
            }
        }
        ret += "}";
        return ret;
    }
};

class Checker {
public:
    ifstream fin1, fin2;
    set<string> s1;
    Checker(string filename1, string filename2) {
        fin1.open(filename1);
        fin2.open(filename2);
        
        if(!fin1 || !fin2) {
            cout << "Input file could not be opened\n";
            exit(0);
        }
    }
    void compare() {
        file1ToSet();
        
        string str_buffer;
        while(!getline(fin2, str_buffer).eof()){
            str_buffer.pop_back();
            if(s1.find(str_buffer) == s1.end()) {
                cout << "failed at " << str_buffer <<  endl;
                return;
            }
        }
    }
    void file1ToSet() {
        string str_buffer;
        while(!getline(fin1, str_buffer).eof()){
            s1.insert(str_buffer);
        }
    }
};

int main (/* int argc, char ** argv */) {
    // if(argc!=4) {
    //     cout << "error : The number of parameters must be 3";
    //     return 0;
    // }
    // string minSupport(argv[1]);
    // string inputFileName(argv[2]);
    // string outputFileName(argv[3]);
    
    /*
    std::vector<Transaction > transactions = {
        {1, 2, 5},
        {2,4},
        {2,3},
        {1, 2, 4},
        {1, 3},
        {2, 3},
        {1, 3},
        {1, 2, 3, 5},
        {1, 2, 3}
    };
    */
    
    std::cout << "Initializing datasets...\n";
    InputReader inputReader("chess.dat");
    std::vector<Transaction> transactions = inputReader.getTransactions();
    std::cout << "Processing Apriori...\n";
    Apriori apriori(transactions, stold("0.2"));
    apriori.process();
    // std::cout << "Output generation...\n";
    // OutputPrinter outputPrinter("out.dat", apriori.getAssociationRules());
    
    /*
    for test
    Checker checker("output5.txt", "outputRsupport5.txt");
    checker.compare();
    */

    
    return 0;
}