#include "FastaMatcher.h"

int main() {
    FastaMatcher mySequences;

    std::string filePath = "C:/Users/kamis/OneDrive/Pulpit/c++/cpp obiektowe/sequence_test.txt";  

    mySequences.readFasta(filePath);

    std::string motif = "ACGTAG";
    int k = 2; 

    for (const auto& sequence : mySequences.getSequences()) {
        std::string seq = sequence.getSeq();

        mySequences.find_matches_with_k_differences(motif, seq, k);
    }

    return 0;
}
