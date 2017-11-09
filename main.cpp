#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <map>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Sequence.h"
#include "SequenceReader.h"

using namespace std;

typedef vector<map<string, double> > Profile;
typedef vector<Sequence> Seqs;
typedef vector<string> Words;
typedef map<string, string> Motifs;

void    printUSAGE    ();
void    printMotifs   (const Motifs& motifs);
void    printProfile  (const Profile& profile, const Words& nwords);
void    printState    (const Motifs& motifs, const Profile& profile, const double entropy,
                       const string& header, const Words& nwords);
void    makeNwords    (const Words& elems, Words* target, string str, int n, int step);
Profile makeProfile   (const Motifs& motifs, const int k, const Words& nwords, int bit);
Profile initiate      (const Seqs& seqs, const int k, Motifs& motifs, const Words& nwords, int bit);
void    expect        (const Seqs& seqs, const int n, const int k, Motifs& motifs, const Profile& profile, int bit);
double  Entropy       (const Profile& profile);
pair<Profile, Motifs>
        run           (const Seqs& seqs, const int k, const Words& nwords, int bit);

int main(int argc, char* argv[]) {
    if (argc < 5 || strcmp(argv[1], "-h") == 0) printUSAGE();

    srand(time(NULL));

    // argument parsing
    int k = 0, step = 100, bit = 1;
    SequenceReader reader;

    for (int i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-k") == 0)      k = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "-i") == 0) reader.open(argv[i + 1]);
        else if (strcmp(argv[i], "-s") == 0) step = atoi(argv[i+1]);
        else if (strcmp(argv[i], "-b") == 0) bit  = atoi(argv[i+1]);
    }
	
	if (!reader.is_open()) {
		cerr << "Input file was not opened." << endl;
		return -1;
	}

    auto seqs = reader.getAll();

    // generate a set of words
    vector<string> nwords, elems = {"A", "C", "G", "T"};
    // element vector정의
    // for (char c = 0, c < 256; ++c) elems.push_back(to_string(c));
    makeNwords(elems, &nwords, "", bit, 0);

    // Input 정보 출력
    cout << "encoding:  " << bit  << endl;
    cout << "nwords:    {";
    for (auto& word : nwords) cout << " " << word;
    cout << " }" << endl;
    cout << "step:      " << step << endl;
    cout << "k:         " << k << endl;
    cout << "# of seqs: " << seqs.size() << endl;

    Profile result_profile;
    Motifs  result_motifs;
    double  result_entropy = -log2(1.0/nwords.size()) * k;
    bool    converge = false;

    // step size를 한 단위로 값이 수렴할 떄 까지 EM 수행
    for (int run_cnt = 0; true; ++run_cnt) {
        if (run_cnt % step == 0 && converge) break;
        else if (run_cnt % step == 0) converge = true;

        auto local_max     = run(seqs, k, nwords, bit);
        auto local_entropy = Entropy(local_max.first);

        // 이전 step 보다 좋은 결과가 나왔을 경우 수렴X
        if (local_entropy < result_entropy) {
            result_entropy = local_entropy;
            result_profile = local_max.first;
            result_motifs  = local_max.second;
            converge       = false;
        }
    }

    printState(result_motifs, result_profile, result_entropy, "[GLOBAL BEST]", nwords);

    return 0;
}

// 사용법을 print하는 함수
void printUSAGE() {
    cout << "MoFi (Motif Finder) v0.0" << endl;
    cout << "usage: ./mofi [-h] [-k #] [-i fasta]" << endl;
    cout << "       -h:\tprint this usage" << endl;
    cout << "       -b:\tencoding bit(s), default 1" << endl;
    cout << "       -s:\tstep run size, default 100" << endl;
    cout << "       -k #\tmotif length" << endl;
    cout << "       -i fasta\tinput fasta file" << endl;
}

// 각 시퀀스에서 찾아진 Motif를 출력하는 함수
void printMotifs(const Motifs& motifs) {
    cout << "# Motifs" << endl;
    for (auto& motif : motifs) {
        cout << motif.first << ": " << motif.second << endl;
    }
    cout << endl;
}

// 전체  Motif의 profile을 출력하는 함수
void printProfile(const Profile& profile, const Words& nwords) {
    for (int i = 0; i < 3 + nwords[0].size(); ++i) cout << " ";
    for (int i = 0; i < profile.size(); ++i) {
        string pos = "[" + to_string(i) + "]";
        cout << setw(10) << pos;
    }
    cout << endl;

    // [+] nbyte 단위로 encoding할 경우 word를 byte로 출력하도록 수정하는 것이 좋음
    for (auto word : nwords) {
        cout << "[" << word << "] ";
        for (auto& prob : profile) cout << setprecision(3) << setw(10) << prob.at(word);
        cout << endl;
    }
    cout << endl;
}

// motif와 그에 해당하는 profile을 차례로 출력
void printState(const Motifs& motifs, const Profile& profile, const double entropy,
                const string& header, const Words& nwords) {
    cout << header << endl;
    printMotifs(motifs);
    printProfile(profile, nwords);
    cout << "Entropy: " << entropy << endl;
}

// element를 기준으로 모든 조합을 만드는 함수
void makeNwords(const Words& elems, Words* target, string str, int n, int step) {
    if (step == n) target->push_back(str);
    else for (auto elem: elems) makeNwords(elems, target, str + elem, n, step + 1);
}

/*
 * 주어진 motif를 기반으로 profile을 만드는 함수
 * profile[pos][word]: motif의 해당 pos위치에서 word가 등장할 확률
 * EM에서 Maximization에 해당
 */
Profile makeProfile(const Motifs& motifs, const int k, const Words& nwords, int bit) {
    Profile profile;
    for (int i = 0; i < k; ++i) {
        map<string, double> prob;
        // profile에 0이 포함되지 않게 하기 위한 pseudo count
        for (auto& word : nwords) prob[word] = 1.0;
        for (auto& motif : motifs) prob[motif.second.substr(i*bit, bit)]++;
        for (auto& item : prob) item.second /= (motifs.size() + nwords.size());

        profile.push_back(prob);
    }

    return profile;
}

/*
 * EM 알고리즘에서 init step에 해당
 * 입력 시퀀스에서 랜덤하게 k-mer를 추출
 */
Profile initiate(const Seqs& seqs, const int k, Motifs& motifs, const Words& nwords, int bit) {
    for (auto& seq : seqs) {
        auto s_pos = rand() % (seq.length() - k * bit);
        motifs[seq.getHeader()] = seq.sub_seq(s_pos, s_pos + k * bit);
    }
    return makeProfile(motifs, k, nwords, bit);
}

/*
 * EM 알고리즘에서 expectation에 해당
 * target sequence에서 나올 수 있는 모든 k-mer에 대해 profile에 기반한 확률 계산
 * 계산한 확률을 기반으로 랜덤하게 k-mer를 추출 (simulated annealing)
 * [+] 확률이 가장 큰 k-mer를 선택하게 수정할 수도 있지만 init step에 영향을 많이 받을 가능성이 있음
 * [+] 가장 확률이 높은 k-mer뿐만 아니라 별로 차이가 안나는 secondary k-mer가 있을 경우 함께 고려가능하다는 장점이 있음
 */
void expect(const Seqs& seqs, const int n, const int k, Motifs& motifs, const Profile& profile, int bit) {
    auto target = seqs.begin() + n;
    // profile biased random select
    vector<double> prob_vec;
    for (int i = 0; i < target->length() - k * bit; ++i) {
        auto kmer = target->sub_seq(i, i + k * bit);

        double prob = 1.0;
        for (size_t pos = 0, step = 0; step < k; pos += bit, ++step)
            prob *= (profile.begin() + step)->at(kmer.substr(pos, bit));
        prob_vec.push_back(prob);
    }

    default_random_engine generator;
    discrete_distribution<int> random_generator(prob_vec.begin(), prob_vec.end());

    auto select = random_generator(generator);

    motifs[target->getHeader()] = target->sub_seq(select, select + k * bit);
}

/*
 * profile의 entropy를 계산
 * entropy가 EM의 heuristic 값 역할
 */
double Entropy(const Profile& profile) {
    double entropy = 0.0;
    for (auto& pos : profile) {
        for (auto& prob : pos) {
            entropy -= prob.second * log2(prob.second);
        }
    }

    return entropy;
}

/*
 * step size만큼 EM을 수행
 * expectation step에서 모든 시퀀스에 대해 expectation을 수행하는 것이 아니라
 * 하나의 시퀀스만 선택해서 expectation을 수행하는 gibbsSampling을 적용
 */
pair<Profile, Motifs> run(const Seqs& seqs, const int k, const Words& nwords, int bit) {
    Motifs  motifs, best_motifs;
    Profile best_profile;
    auto    profile      = initiate(seqs, k, motifs, nwords, bit);
    double  best_entropy = Entropy(profile);

    default_random_engine         generator;
    uniform_int_distribution<int> gibbs_sampler(0, (int)seqs.size() - 1);


//    printState(motifs, profile, best_entropy, "[Init]", nwords);

    int step = 0;
    bool run = true;
    while (run) {
        double  local_best = -log2(1.0/nwords.size()) * k;
        Motifs  local_best_motifs;
        Profile local_best_profile;
        for (int i = 0; i < seqs.size(); ++i) {
            int sampling = gibbs_sampler(generator);

            motifs.erase((seqs.begin() + sampling)->getHeader());
            profile = makeProfile(motifs, k, nwords, bit);
            expect(seqs, sampling, k, motifs, profile, bit);

            double entropy = Entropy(profile);
            if (entropy < local_best) {
                local_best = entropy;
                local_best_motifs = motifs;
                local_best_profile = profile;
            }
        }
//        printState(local_best_motifs, local_best_profile, local_best, "[Step" + to_string(step) + "]", nwords);
//        cout << "Local max: " << best_entropy << endl;
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

    printState(best_motifs, best_profile, best_entropy, "[LOCAL MAX]", nwords);

    return make_pair(best_profile, best_motifs);
}