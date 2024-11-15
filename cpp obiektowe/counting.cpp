#include "FastaMatcher.h"

int main() {
    FastaMatcher matcher;

    std::string fastaFilePath = "C:/Users/kamis/OneDrive/Pulpit/c++/cpp obiektowe/rosalind_dna(2).txt";

    matcher.readFasta(fastaFilePath);

    std::vector<FastaMatcher> sequences = matcher.getSequences();
    
    for (const auto& sequence : sequences) {
        std::string seq = sequence.getSeq(); 
        matcher.count_nucleotides(seq);
    }

    return 0;
}
