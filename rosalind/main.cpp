#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include "rosalind.h"

using namespace std;

int main() {
    int choice;
    cout << "Wybierz zadanie (1-10):" << endl;
    cout << "1. DNACounter - liczenie nukleotydów" << endl;
    cout << "2. Transkrypcja - DNA -> RNA" << endl;
    cout << "3. Translacja - RNA -> białko" << endl;
    cout << "4. Reverse complement" << endl;
    cout << "5. SpliceAndTranslate - usunięcie intronów i translacja" << endl;
    cout << "6. HammingDistanceCalculator - obliczanie odległości Hamminga" << endl;
    cout << "7. SubstringLocator - wyszukiwanie motywu w sekwencji" << endl;
    cout << "8. SubsequenceFinder - wyszukiwanie indeksów subsekwencji" << endl;
    cout << "9. LongestCommonSubsequenceFinder - najdłuższa wspólna subsekwencja" << endl;
    cout << "10. GeneralEditDistanceCalculator - obliczanie odległości edycji" << endl;
    cout << "Podaj numer zadania: ";
    cin >> choice;
    cin.ignore();

    switch(choice) {
        case 1: {
            cout << "Podaj sekwencję:" << endl;
            string seq = readSequence();
            DNACounter counter(seq);
            int countA, countC, countG, countT;
            counter.countNucleotides(countA, countC, countG, countT);
            cout << "A: " << countA << ", C: " << countC 
                 << ", G: " << countG << ", T: " << countT << endl;
            break;
        }
        case 2: {
            cout << "Podaj sekwencję:" << endl;
            string seq = readSequence();
            SequenceAnalyzer analyzer(seq);
            cout << "RNA: " << analyzer.transcribe() << endl;
            break;
        }
        case 3: {
            cout << "Podaj sekwencję RNA:" << endl;
            string rna = readSequence();
            SequenceAnalyzer analyzer(rna);
            cout << "Białko: " << analyzer.translate() << endl;
            break;
        }
        case 4: {
            cout << "Podaj sekwencję:" << endl;
            string seq = readSequence();
            SequenceAnalyzer analyzer(seq);
            cout << "Reverse complement: " << analyzer.getReverseComplement() << endl;
            break;
        }
        case 5: {
            cout << "Podaj sekwencje w formacie FASTA:" << endl;
            vector<pair<string, string>> fasta = readFasta();
            if(fasta.empty()){
                break;
            }
            string mainSeq = fasta[0].second;
            vector<GeneticSequence> introns;
            for (size_t i = 1; i < fasta.size(); i++){
                introns.push_back(GeneticSequence(fasta[i].second));
            }
            SpliceAndTranslate st(mainSeq, introns);
            cout << "Protein string: " << st.spliceTranslate() << endl;
            break;
        }        
        case 6: {
            cout << "Podaj pierwszą sekwencję:" << endl;
            string s = readSequence();
            cout << "Podaj drugą sekwencję:" << endl;
            string t = readSequence();
            if(s.size() != t.size()){
                cout << "Sekwencje muszą mieć równą długość!" << endl;
                break;
            }
            HammingDistanceCalculator hdc{ GeneticSequence(s), GeneticSequence(t) };
            int distance = hdc.calculate();
            if(distance < 0)
                cout << "Błąd: sekwencje różnej długości." << endl;
            else
                cout << "Odległość Hamminga: " << distance << endl;
            break;
        }
        case 7: {
            cout << "Podaj sekwencję::" << endl;
            string mainSeq = readSequence();
            cout << "Podaj motyw:" << endl;
            string motif = readSequence();
            SubstringLocator locator{ GeneticSequence(mainSeq), GeneticSequence(motif) };
            vector<int> locations = locator.findLocations();
            if(locations.empty())
                cout << "Motyw nie został znaleziony." << endl;
            else {
                cout << "Lokalizacje motywu:" << endl;
                for (int pos : locations)
                    cout << pos << " ";
                cout << endl;
            }
            break;
        }
        case 8: {
            cout << "Podaj sekwencję:" << endl;
            string mainSeq = readSequence();
            cout << "Podaj subsekwencję:" << endl;
            string subseq = readSequence();
            SubsequenceFinder finder{ GeneticSequence(mainSeq), GeneticSequence(subseq) };
            vector<int> indices = finder.findIndices();
            if(indices.empty())
                cout << "Subsekwencja nie została znaleziona w głównej sekwencji." << endl;
            else {
                cout << "Indeksy występowania kolejnych znaków subsekwencji:" << endl;
                for (int idx : indices)
                    cout << idx << " ";
                cout << endl;
            }
            break;
        }
        case 9: {
            cout << "Podaj sekwencje w formacie FASTA:" << endl;
            vector<pair<string, string>> fasta = readFasta();
            if(fasta.size() < 2) {
                break;
            }
            LongestCommonSubsequenceFinder lcsFinder{ GeneticSequence(fasta[0].second),
                                                       GeneticSequence(fasta[1].second) };
            string lcs = lcsFinder.findLongestCommonSubsequence();
            if(lcs.empty())
                cout << "Brak wspólnej subsekwencji." << endl;
            else {
                cout << "Najdłuższa wspólna subsekwencja:" << endl;
                cout << lcs << endl;
            }
            break;
        }
        case 10: {
            cout << "Podaj sekwencje w formacie FASTA:" << endl;
            vector<pair<string, string>> fasta = readFasta();
            if(fasta.size() < 2) {
                break;
            }
            GeneralEditDistanceCalculator ged{ GeneticSequence(fasta[0].second),
                                                 GeneticSequence(fasta[1].second) };
            cout << "Edit distance: " << ged.calculate() << endl;
            break;
        }
        default:
            cout << "Nieprawidłowy wybór!" << endl;
    }
    return 0;
}
