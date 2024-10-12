#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

class fastaseq {
public:
    string id;   
    string seq;  

    vector<fastaseq> sequences;

    void readFasta(const string& filePath) {
        ifstream file(filePath);
        
        if (!file.is_open()) {
            cerr << "Error: Could not open the file: " << filePath << endl;
            return;
        }

        string line;
        fastaseq currentSeq;  
        while (getline(file, line)) {
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

};

int main() {
    fastaseq mySequences; 

    string filePath = "C:/Users/kamis/OneDrive/Pulpit/c++/cpp obiektowe/sequence_test.txt";  

    mySequences.readFasta(filePath);


    return 0;
}