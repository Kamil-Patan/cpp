#include "FastaMatcher.h"

void FastaMatcher::readFasta(const std::string& filePath) {
    std::ifstream file(filePath);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file: " << filePath << std::endl;
        return;
    }

    std::string line;
    FastaMatcher currentSeq;  
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!currentSeq.id.empty()) {
                sequences.push_back(currentSeq);  
                currentSeq = FastaMatcher();          
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

int FastaMatcher::count_differences(const std::string &motif, const std::string &seq) {
    int differences = 0;
    for (size_t i = 0; i < motif.length(); ++i) {
        if (motif[i] != seq[i]) {
            differences++;
        }
    }
    return differences;
}

void FastaMatcher::find_matches_with_k_differences(const std::string &motif, const std::string &seq, int k) {
    size_t len_motif = motif.length();
    std::vector<std::pair<int, int>> found_matches;

    for (size_t i = 0; i <= seq.length() - len_motif; ++i) {
        std::string fragment = seq.substr(i, len_motif);  
        
        int differences = count_differences(motif, fragment);
        
        if (differences <= k) {
            found_matches.emplace_back(i + 1, i + len_motif); 
        }
    }

    if (!found_matches.empty()) {
        for (const auto &match : found_matches) {
            int start = match.first;
            int end = match.second;

            for (int j = end - k; j <= end; ++j) {
                if (j >= start) { 
                    std::cout << start << " " << j << std::endl;
                }
            }
        }
    }
}
