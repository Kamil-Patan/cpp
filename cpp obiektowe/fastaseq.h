// fastaseq.h
#ifndef FASTASEQ_H
#define FASTASEQ_H

#include <string>
#include <vector>

class fastaseq {
public:
    std::string id;   
    std::string seq;  

    std::vector<fastaseq> sequences;

    void readFasta(const std::string& filePath);
};

#endif // FASTASEQ_H
