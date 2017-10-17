//
// Created by hjgwak(Kelvin) on 07/10/2016.
//

#include "Sequence.h"

// Constructor
Sequence::Sequence() {
    this->header = "";
    this->sequence = "";
}

Sequence::Sequence(const std::string &header, const std::string &sequence) {
    this->header = header;
    this->sequence = sequence;
}

Sequence::Sequence(const std::string& header, const char* sequence) {
	this->header = header;
	this->sequence = sequence;
}

Sequence::Sequence(const std::string& header, const char ch) {
	this->header = header;
	this->sequence = ch;
}

Sequence::Sequence(const Sequence &seq) {
    this->header = seq.header;
    this->sequence = seq.sequence;
}

Sequence &Sequence::operator=(const std::string &sequence) {
    this->sequence = sequence;

    return *this;
}

Sequence &Sequence::operator=(const char *sequence) {
    this->sequence = sequence;

    return *this;
}

Sequence &Sequence::operator=(const char c) {
    this->sequence = c;

    return *this;
}

Sequence &Sequence::operator=(const Sequence &sequence) {
    this->header = sequence.header;
    this->sequence = sequence.sequence;

    return *this;
}

// Setters

void Sequence::setHeader(const std::string &header) {
    this->header = header;
}

void Sequence::append(const std::string &seq) {
    this->sequence += seq;
}

void Sequence::append(const Sequence &seq) {
    this->sequence += seq.sequence;
}

void Sequence::append(const char *seq) {
    this->sequence += seq;
}

void Sequence::append(const char c) {
    this->sequence += c;
}

// Getters

std::string Sequence::getHeader() const {
    return this->header;
}

std::string Sequence::seq() const {
    return this->sequence;
}

std::string Sequence::sub_seq(const int start, const int end) const {
    int real_start = (start < 0) ? (int)this->length() + start : start;
    int real_end = (end <= 0) ? (int)this->length() + end : end;

    if (real_start < 0 || real_end < 0) {
        return "\0";
    }

    int len = real_end - real_start;
    if (len < 0) return "\0";

    return this->sequence.substr((size_t)real_start, (size_t)len);
}

size_t Sequence::length() const {
    return this->sequence.size();
}

char null = '\0';

char &Sequence::operator[](const int pos) {
    int real_pos = (pos < 0) ? (int)this->length() + pos : pos;

    if (real_pos > this->length() || real_pos < 0) return null;
    else return this->sequence[real_pos];
}

const char &Sequence::operator[](const int pos) const {
    int real_pos = (pos < 0) ? (int)this->length() + pos : pos;

    if (real_pos > this->length() || real_pos < 0) return null;
    else return this->sequence[real_pos];
}

// Operators

Sequence &Sequence::operator+=(const std::string &sequence) {
    this->append(sequence);
    return *this;
}

Sequence &Sequence::operator+=(const char *sequence) {
    this->append(sequence);
    return *this;
}

Sequence &Sequence::operator+=(const char c) {
    this->append(c);
    return *this;
}

Sequence &Sequence::operator+=(const Sequence &sequence) {
    this->append(sequence);
    return *this;
}

Sequence Sequence::operator+ (const Sequence& rhs) {
    Sequence new_seq(*this);
    new_seq += rhs;

    return new_seq;
}

Sequence Sequence::operator+ (const std::string& rhs) {
    Sequence new_seq(*this);
    new_seq += rhs;

    return new_seq;
}

Sequence Sequence::operator+ (const char* rhs) {
    Sequence new_seq(*this);
    new_seq += rhs;

    return new_seq;
}

Sequence Sequence::operator+ (const char rhs) {
    Sequence new_seq(*this);
    new_seq += rhs;

    return new_seq;
}

Sequence operator+ (const std::string& lhs, const Sequence& rhs) {
    Sequence new_seq(rhs.getHeader(), lhs);
    new_seq += rhs;

    return new_seq;
}

Sequence operator+ (const char* lhs, const Sequence& rhs) {
    Sequence new_seq(rhs.getHeader(), lhs);
    new_seq += rhs;

    return new_seq;
}

Sequence operator+ (const char lhs, const Sequence& rhs) {
    Sequence new_seq(rhs.getHeader(), lhs);
    new_seq += rhs;

    return new_seq;
}
