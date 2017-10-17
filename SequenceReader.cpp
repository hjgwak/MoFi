#include "SequenceReader.h"

using namespace std;

//constructors && distructors
SequenceReader::SequenceReader() {

}

SequenceReader::SequenceReader(const char* file_name) {
	this->open(file_name);
}

SequenceReader::SequenceReader(const string& file_name) {
	this->open(file_name);
}

SequenceReader::~SequenceReader() {
	this->close();
}

//setters
void SequenceReader::open(const char* file_name) {
	if (!this->file.is_open()) {
		this->file.open(file_name, ifstream::in);
	}
}

void SequenceReader::open(const string& file_name) {
	if (!this->file.is_open()) {
		this->file.open(file_name, ifstream::in);
	}
}

void SequenceReader::close() {
	if (this->file.is_open()) {
		this->file.close();
	}
}

//getters
Sequence SequenceReader::getNext() {
	if (!this->hasNext()) {
		return Sequence();
	}

	if (this->iqueue.empty()) {
		this->readNext();
	}

	Sequence seq = this->iqueue.front();
	this->iqueue.pop();

	return seq;
}

vector<Sequence> SequenceReader::getAll() {
	vector<Sequence> seqs;
	this->readAll();

	while (!this->iqueue.empty()) {
		seqs.push_back(this->iqueue.front());
		this->iqueue.pop();
	}

	return seqs;
}

bool SequenceReader::hasNext() const {
	return !this->iqueue.empty() || this->file.good();
}

bool SequenceReader::is_open() const {
	return this->file.is_open();
}

//actions
bool SequenceReader::readNext() {
	if (this->file.good()) {
		string line, header, seq;
		streampos prev;

		seq = "";

		getline(this->file, header);
		while (this->file.good()) {
			getline(this->file, line);
			if (line[0] == '>') {
				this->file.seekg(prev);
				break;
			}
			else {
				seq += line;
				prev = this->file.tellg();
			}
		}

		this->iqueue.push(Sequence(header, seq));

		return true;
	}

	return false;
}

bool SequenceReader::readAll() {
	while (this->file.good()) {
		this->readNext();
	}

	return false;
}