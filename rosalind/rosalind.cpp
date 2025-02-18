#include "rosalind.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <vector>
#include <limits>
#include <unordered_map>
using namespace std;

string readSequence() {
    string seq;
    getline(cin, seq);
    while(seq.empty() && !cin.eof()){
        getline(cin, seq);
    }
    return seq;
}

vector<pair<string, string>> readFasta() {
    vector<pair<string, string>> fasta;
    string line, id, seq;
    while(getline(cin, line)) {
        if(line.empty()) break;
        if(line[0] == '>') {
            if(!id.empty()){
                fasta.push_back({id, seq});
            }
            id = line.substr(1);
            seq = "";
        } else {
            seq += line;
        }
    }
    if(!id.empty()){
        fasta.push_back({id, seq});
    }
    return fasta;
}

GeneticSequence::GeneticSequence(const string &seq) {
    setSequence(seq);
}

const string &GeneticSequence::getSequence() const {
    return sequence;
}

void GeneticSequence::setSequence(const string &seq) {
    sequence = seq;
}


SequenceAnalyzer::SequenceAnalyzer(const string &seq) : GeneticSequence(seq) {}

string SequenceAnalyzer::transcribe() const {
    string rna = getSequence();
    for (char &c : rna) {
        if(c == 'T') c = 'U';
    }
    return rna;
}

string SequenceAnalyzer::translate() const {
    static const unordered_map<string, char> codonTable = {
        {"UUU", 'F'}, {"UUC", 'F'},
        {"UUA", 'L'}, {"UUG", 'L'},
        {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
        {"UAU", 'Y'}, {"UAC", 'Y'},
        {"UAA", ' '}, {"UAG", ' '}, {"UGA", ' '},
        {"UGU", 'C'}, {"UGC", 'C'},
        {"UGG", 'W'},
        {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
        {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAU", 'H'}, {"CAC", 'H'},
        {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'},
        {"AUG", 'M'},
        {"ACU", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAU", 'N'}, {"AAC", 'N'},
        {"AAA", 'K'}, {"AAG", 'K'},
        {"AGU", 'S'}, {"AGC", 'S'},
        {"AGA", 'R'}, {"AGG", 'R'},
        {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'}, {"GUG", 'V'},
        {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAU", 'D'}, {"GAC", 'D'},
        {"GAA", 'E'}, {"GAG", 'E'},
        {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
    
    string rna = transcribe();
    string protein;
    for (size_t i = 0; i + 2 < rna.size(); i += 3) {
        string codon = rna.substr(i, 3);
        if (codon == "UAA" || codon == "UAG" || codon == "UGA")
            break;
        auto it = codonTable.find(codon);
        if (it != codonTable.end())
            protein.push_back(it->second);
    }
    return protein;
}

string SequenceAnalyzer::getReverseComplement() const {
    string revComp;
    string seq = getSequence();
    for(auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char nucleotide = *it;
        char comp;
        switch(nucleotide) {
            case 'A': comp = 'T'; break;
            case 'T': comp = 'A'; break;
            case 'C': comp = 'G'; break;
            case 'G': comp = 'C'; break;
            default: comp = nucleotide; break;
        }
        revComp.push_back(comp);
    }
    return revComp;
}

SpliceAndTranslate::SpliceAndTranslate(const string &mainSeq, const vector<GeneticSequence> &introns)
    : SequenceAnalyzer(mainSeq), introns(introns) {}

string SpliceAndTranslate::spliceTranslate() const {
    string exon = getSequence();
    for (const auto &intron : introns) {
        string intronSeq = intron.getSequence();
        size_t pos = exon.find(intronSeq);
        while(pos != string::npos) {
            exon.erase(pos, intronSeq.length());
            pos = exon.find(intronSeq);
        }
    }
    SequenceAnalyzer analyzer(exon);
    return analyzer.translate();
}

PairwiseSequenceOperation::PairwiseSequenceOperation(const GeneticSequence &s1, const GeneticSequence &s2)
    : seq1(s1), seq2(s2) {}

const string &PairwiseSequenceOperation::getSequence1() const {
    return seq1.getSequence();
}

const string &PairwiseSequenceOperation::getSequence2() const {
    return seq2.getSequence();
}

DNACounter::DNACounter(const string &seq) : GeneticSequence(seq) {}

void DNACounter::countNucleotides(int &countA, int &countC, int &countG, int &countT) const {
    countA = countC = countG = countT = 0;
    for (char nucleotide : getSequence()) {
        switch(nucleotide) {
            case 'A': countA++; break;
            case 'C': countC++; break;
            case 'G': countG++; break;
            case 'T': countT++; break;
            default: break;
        }
    }
}

HammingDistanceCalculator::HammingDistanceCalculator(const GeneticSequence &s, const GeneticSequence &t)
    : PairwiseSequenceOperation(s, t) {}

int HammingDistanceCalculator::calculate() const {
    const string &s = getSequence1();
    const string &t = getSequence2();
    if(s.size() != t.size())
        return -1;
    int distance = 0;
    for (size_t i = 0; i < s.size(); i++) {
        if(s[i] != t[i])
            distance++;
    }
    return distance;
}

SubstringLocator::SubstringLocator(const GeneticSequence &mainSeq, const GeneticSequence &pattern)
    : PairwiseSequenceOperation(mainSeq, pattern) {}

vector<int> SubstringLocator::findLocations() const {
    vector<int> locations;
    const string &s = getSequence1();
    const string &pattern = getSequence2();
    size_t pos = 0;
    while((pos = s.find(pattern, pos)) != string::npos) {
        locations.push_back(pos + 1);
        pos++;
    }
    return locations;
}

SubsequenceFinder::SubsequenceFinder(const GeneticSequence &mainSeq, const GeneticSequence &subseq)
    : PairwiseSequenceOperation(mainSeq, subseq) {}

vector<int> SubsequenceFinder::findIndices() const {
    vector<int> indices;
    const string &s = getSequence1();
    const string &subseq = getSequence2();
    size_t pos = 0;
    for (char c : subseq) {
        bool found = false;
        while (pos < s.size()) {
            if(s[pos] == c) {
                indices.push_back(pos + 1);
                pos++;
                found = true;
                break;
            }
            pos++;
        }
        if (!found) {
            indices.clear();
            break;
        }
    }
    return indices;
}

LongestCommonSubsequenceFinder::LongestCommonSubsequenceFinder(const GeneticSequence &s, const GeneticSequence &t)
    : PairwiseSequenceOperation(s, t) {}

string LongestCommonSubsequenceFinder::findLongestCommonSubsequence() const {
    const string &sStr = getSequence1();
    const string &tStr = getSequence2();
    int m = sStr.size(), n = tStr.size();
    vector<vector<int>> dp(m+1, vector<int>(n+1, 0));
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if(sStr[i-1] == tStr[j-1])
                dp[i][j] = dp[i-1][j-1] + 1;
            else
                dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
        }
    }
    int i = m, j = n;
    string lcs;
    while(i > 0 && j > 0) {
        if(sStr[i-1] == tStr[j-1]) {
            lcs.push_back(sStr[i-1]);
            i--; j--;
        } else if(dp[i-1][j] >= dp[i][j-1]) {
            i--;
        } else {
            j--;
        }
    }
    reverse(lcs.begin(), lcs.end());
    return lcs;
}

GeneralEditDistanceCalculator::GeneralEditDistanceCalculator(const GeneticSequence &s, const GeneticSequence &t)
    : PairwiseSequenceOperation(s, t) {}

int GeneralEditDistanceCalculator::calculate() const {
    const string &sStr = getSequence1();
    const string &tStr = getSequence2();
    int m = sStr.size(), n = tStr.size();
    vector<vector<int>> dp(m+1, vector<int>(n+1, 0));
    for (int i = 0; i <= m; i++)
        dp[i][0] = i;
    for (int j = 0; j <= n; j++)
        dp[0][j] = j;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int cost = (sStr[i-1] == tStr[j-1]) ? 0 : 1;
            dp[i][j] = min({ dp[i-1][j] + 1,
                             dp[i][j-1] + 1,
                             dp[i-1][j-1] + cost });
        }
    }
    return dp[m][n];
}
