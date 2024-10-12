#include <iostream>
#include <string>
#include <vector>
#include <fstream>

class FastaReader {
private:
    std::string dnaString;
    std::string motif;

public:
    FastaReader(std::string filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                if (line[0] != '>') {
                    dnaString += line;
                }
            }
            file.close();
        }
    }

    void setMotif(std::string m) {
        motif = m;
    }

    std::vector<std::pair<int, int>> findSubstrings(int k) {
        std::vector<std::pair<int, int>> substrings;
        
        return substrings;
    }
};

int main() {
    FastaReader reader("C:/Users/kamis/OneDrive/Pulpit/c++/cpp obiektowe/sequence_test.txt");
    reader.setMotif("ACGTAG");
    int k = 2;
    std::vector<std::pair<int, int>> result = reader.findSubstrings(k);

    for (const auto& substr : result) {
        std::cout << substr.first << " " << substr.second << std::endl;
    }

    return 0;
}