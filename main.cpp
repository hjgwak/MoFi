#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Sequence.h"
#include "SequenceReader.h"

using namespace std;

typedef vector<map<char, double> > Profile;
typedef vector<Sequence> Seqs;
typedef map<string, string> Motifs;

void printUSAGE();
void printMotifs(const Motifs& motifs);
void printProfile(const Profile& profile);
void printState(const Motifs& motifs, const Profile& profile, const double entropy, const string& header);
Profile makeProfile(const Motifs& motifs, const int k);
Profile initiate(const Seqs& seqs, const int k, Motifs& motifs);
void expect(const Seqs& seqs, const int n, const int k, Motifs& motifs, const Profile& profile);
double Entropy(const Profile& profile);

int main(int argc, char* argv[]) {
    if (argc < 5 || strcmp(argv[1], "-h") == 0) printUSAGE();

    srand(time(NULL));

    int k = 0;
    SequenceReader reader;

    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-k") == 0)      k = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "-i") == 0) reader.open(argv[i + 1]);
    }
	
	if (!reader.is_open()) {
		cerr << "Input file was not opened." << endl;
		system("pause");
		return -1;
	}

    auto seqs = reader.getAll();

    cout << "k:         " << k << endl;
    cout << "# of seqs: " << seqs.size() << endl;

    Motifs  motifs, best_motifs;
    Profile best_profile;
    auto    profile      = initiate(seqs, k, motifs);
    double  best_entropy = Entropy(profile);

    printState(motifs, profile, best_entropy, "[Init]");

    int step = 0;
    bool run = true;
    while (run) {
        double  local_best = 3.0 * k;
        Motifs  local_best_motifs;
        Profile local_best_profile;
        for (int i = 0; i < seqs.size(); ++i) {
            motifs.erase((seqs.begin() + i)->getHeader());
            profile = makeProfile(motifs, k);
            expect(seqs, i, k, motifs, profile);

            double entropy = Entropy(profile);
            if (entropy < local_best) {
                local_best = entropy;
                local_best_motifs = motifs;
                local_best_profile = profile;
            }
            /*
            cout << "[Step" << i << "]" << endl;
            printMotifs(motifs);
            printProfile(profile);
             */
        }
        printState(local_best_motifs, local_best_profile, local_best, "[Step" + to_string(step) + "]");
        cout << "Global Best: " << best_entropy << endl;
        step++;

        if (local_best < best_entropy) {
            best_entropy = local_best;
            best_motifs  = local_best_motifs;
            best_profile = local_best_profile;
        }
        else {
            run          = false;
        }
    }

    printState(best_motifs, best_profile, best_entropy, "[BEST]");

    return 0;
}

void printUSAGE() {
    cout << "MoFi (Motif Finder) v0.0" << endl;
    cout << "usage: ./mofi [-h] [-k #] [-i fasta]" << endl;
    cout << "       -h:\tprint this usage" << endl;
    cout << "       -k #\tmotif length" << endl;
    cout << "       -i fasta\tinput fasta file" << endl;
}

void printMotifs(const Motifs& motifs) {
    cout << "# Motifs" << endl;
    for (auto& motif : motifs) {
        cout << motif.first << ": " << motif.second << endl;
    }
    cout << endl;
}

void printProfile(const Profile& profile) {
    cout << "    ";
    for (int i = 0; i < profile.size(); ++i) {
        string pos = "[" + to_string(i) + "]";
        cout << setw(10) << pos;
    }
    cout << endl;

    for (char ch : {'A', 'C', 'G', 'T'}) {
        cout << "[" << ch << "] ";
        for (auto& prob : profile) cout << setprecision(3) << setw(10) << prob.at(ch);
        cout << endl;
    }
    cout << endl;
}

void printState(const Motifs& motifs, const Profile& profile, const double entropy, const string& header) {
    cout << header << endl;
    printMotifs(motifs);
    printProfile(profile);
    cout << "Entropy: " << entropy << endl;
}

Profile makeProfile(const Motifs& motifs, const int k) {
    Profile profile;
    for (int i = 0; i < k; ++i) {
        /* epsilon */
        /*
        map<char, double> prob = {{'A', 0.0}, {'C', 0.0}, {'G', 0.0}, {'T', 0.0}};
        for (auto& motif : motifs) prob[motif.second[i]]++;
        for (auto& item  : prob)
            item.second = (item.second == 0) ? numeric_limits<double>::min() : item.second / motifs.size();
         */
        /* pseudo count */

        map<char, double> prob = {{'A', 1.0}, {'C', 1.0}, {'G', 1.0}, {'T', 1.0}};
        for (auto& motif : motifs) prob[motif.second[i]]++;
        for (auto& item : prob) item.second /= (motifs.size()+4);

        profile.push_back(prob);
    }

    return profile;
}

Profile initiate(const Seqs& seqs, const int k, Motifs& motifs) {
    for (auto& seq : seqs) {
        auto s_pos = rand() % (seq.length() - k);
        motifs[seq.getHeader()] = seq.sub_seq(s_pos, s_pos + k);
    }

    return makeProfile(motifs, k);
}

void expect(const Seqs& seqs, const int n, const int k, Motifs& motifs, const Profile& profile) {
    auto target = seqs.begin() + n;
    double max_prob = -1.0;
    string max_kmer = "";

    for (int i = 0; i < target->length() - k; ++i) {
        auto kmer = target->sub_seq(i, i + k);

        double prob = 1.0;
        for (size_t pos = 0; pos < k; ++pos) prob *= (profile.begin() + pos)->at(kmer[pos]);

        if (prob > max_prob) {
            max_prob = prob;
            max_kmer = kmer;
        }
    }

    motifs[target->getHeader()] = max_kmer;
}

double Entropy(const Profile& profile) {
    double entropy = 0.0;
    for (auto& pos : profile) {
        for (auto& prob : pos) {
            entropy -= prob.second * log2(prob.second);
        }
    }

    return entropy;
}