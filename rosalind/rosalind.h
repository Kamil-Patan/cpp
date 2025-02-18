#ifndef ROSALIND_H
#define ROSALIND_H

#include <string>
#include <vector>
#include <utility>
using namespace std;


class GeneticSequence {
protected:
    string sequence;
public:
    GeneticSequence(const string &seq);
    const string &getSequence() const;
    void setSequence(const string &seq);
};

string readSequence();
vector<pair<string, string>> readFasta();

class SequenceAnalyzer : public GeneticSequence {
public:
    SequenceAnalyzer(const string &seq);
    string transcribe() const;
    virtual string translate() const;
    string getReverseComplement() const;
};

class SpliceAndTranslate : public SequenceAnalyzer {
private:
    vector<GeneticSequence> introns;
public:
    SpliceAndTranslate(const string &mainSeq, const vector<GeneticSequence> &introns);
    string spliceTranslate() const;
};

class PairwiseSequenceOperation {
protected:
    GeneticSequence seq1;
    GeneticSequence seq2;
public:
    PairwiseSequenceOperation(const GeneticSequence &s1, const GeneticSequence &s2);
    const string &getSequence1() const;
    const string &getSequence2() const;
};

class DNACounter : public GeneticSequence {
public:
    DNACounter(const string &seq);
    void countNucleotides(int &countA, int &countC, int &countG, int &countT) const;
};

class HammingDistanceCalculator : public PairwiseSequenceOperation {
public:
    HammingDistanceCalculator(const GeneticSequence &s, const GeneticSequence &t);
    int calculate() const;
};

class SubstringLocator : public PairwiseSequenceOperation {
public:
    SubstringLocator(const GeneticSequence &mainSeq, const GeneticSequence &pattern);
    vector<int> findLocations() const;
};

class SubsequenceFinder : public PairwiseSequenceOperation {
public:
    SubsequenceFinder(const GeneticSequence &mainSeq, const GeneticSequence &subseq);
    vector<int> findIndices() const;
};

class LongestCommonSubsequenceFinder : public PairwiseSequenceOperation {
public:
    LongestCommonSubsequenceFinder(const GeneticSequence &s, const GeneticSequence &t);
    string findLongestCommonSubsequence() const;
};

class GeneralEditDistanceCalculator : public PairwiseSequenceOperation {
public:
    GeneralEditDistanceCalculator(const GeneticSequence &s, const GeneticSequence &t);
    int calculate() const;
};

#endif // ROSALIND_H
