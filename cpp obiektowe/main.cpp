// main.cpp
#include "fastaseq.h"
#include <iostream>
#include <string>

int main() {
    fastaseq mySequences; 

    std::string filePath = "C:/Users/kamis/OneDrive/Pulpit/c++/cpp obiektowe/sequence_test.txt";  

    mySequences.readFasta(filePath);

    return 0;
}
