// fastaseq.cpp
#include "fastaseq.h"
#include <fstream>
#include <iostream>

void fastaseq::readFasta(const std::string& filePath) {
    std::ifstream file(filePath);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file: " << filePath << std::endl;
        return;
    }

    std::string line;
    fastaseq currentSeq;  
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!currentSeq.id.empty()) {
                sequences.push_back(currentSeq);  
                currentSeq = fastaseq();          
            }
            currentSeq.id = line.substr(1);
        } else {
            currentSeq.seq += line;  
        }
    }

    if (!currentSeq.id.empty()) {
        sequences.push_back(currentSeq);
    }

    file.close();
}
