#ifndef FASTAMATCHER_H
#define FASTAMATCHER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

class FastaMatcher {
public:
    std::string id;   
    std::string seq;  

    std::vector<FastaMatcher> sequences;

    void readFasta(const std::string& filePath);

    int count_differences(const std::string &motif, const std::string &seq);

    void find_matches_with_k_differences(const std::string &motif, const std::string &seq, int k);
};

#endif
