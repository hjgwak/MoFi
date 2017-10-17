#ifndef __SEQUENCEREADER_H__
#define __SEQUENCEREADER_H__

#include "Sequence.h"
#include <fstream>
#include <queue>
#include <vector>

class SequenceReader {
public :
	//constructors
	SequenceReader();
	SequenceReader(const char* file_name);
	SequenceReader(const std::string& file_name);
	~SequenceReader();

	//setters
	void open(const char* file_name);
	void open(const std::string& file_name);
	void close();

	//getters
	Sequence getNext();
	std::vector<Sequence> getAll();
	bool hasNext() const;
	bool is_open() const;

private :
	//actions
	bool readNext();
	bool readAll();

	std::ifstream file;
	std::queue<Sequence> iqueue;
};

#endif // !__SEQUENCEREADER_H__
