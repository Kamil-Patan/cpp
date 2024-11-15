#ifndef FASTAMATCHER_H
#define FASTAMATCHER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

class FastaMatcher {
private:
    std::string id;
    std::string seq;
    std::vector<FastaMatcher> sequences;

public:
    std::string getId() const { return id; }

    void setId(const std::string& newId) { id = newId; }

    std::string getSeq() const { return seq; }

    void setSeq(const std::string& newSeq) { seq = newSeq; }

    std::vector<FastaMatcher> getSequences() const { return sequences; }

    void setSequences(const std::vector<FastaMatcher>& newSequences) { sequences = newSequences; }

    void readFasta(const std::string& filePath);

    int count_differences(const std::string &motif, const std::string &seq);

    void find_matches_with_k_differences(const std::string &motif, const std::string &seq, int k);

    void count_nucleotides(const std::string &seq);
};

#endif
