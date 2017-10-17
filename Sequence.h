//
// Created by hjgwak(Kelvin) on 07/10/2016.
//

#ifndef SCL_SEQUENCE_H
#define SCL_SEQUENCE_H

#include <string>

class Sequence {
public:
    // Constructors
    Sequence();
    Sequence(const std::string& header, const std::string& sequence);
	Sequence(const std::string& header, const char* sequence);
	Sequence(const std::string& header, const char ch);
    Sequence(const Sequence& seq);

    virtual Sequence& operator=(const std::string& sequence);
    virtual Sequence& operator=(const char* sequence);
    virtual Sequence& operator=(const char c);
    virtual Sequence& operator=(const Sequence& sequence);

    // Setter
    void setHeader(const std::string& header);
    virtual void append(const std::string& seq);
    virtual void append(const Sequence& seq);
    virtual void append(const char* seq);
    virtual void append(const char c);

    // Getter
    std::string getHeader() const ;
    std::string seq() const ;
    std::string sub_seq(const int start, const int end) const;
    size_t length() const;

    char& operator[] (const int pos);
    const char& operator[] (const int pos) const;

    // operator
    virtual Sequence& operator+=(const std::string& sequence);
    virtual Sequence& operator+=(const char* sequence);
    virtual Sequence& operator+=(const char c);
    virtual Sequence& operator+=(const Sequence& sequence);

    Sequence operator+(const Sequence& rhs);
    Sequence operator+(const std::string& rhs);
    Sequence operator+(const char* rhs);
    Sequence operator+(const char rhs);
protected:
    std::string header;
    std::string sequence;
};

// Another operators
Sequence operator+ (const std::string& lhs, const Sequence& rhs);
Sequence operator+ (const char* lhs, const Sequence& rhs);
Sequence operator+ (const char lhs, const Sequence& rhs);

#endif //SCL_SEQUENCE_H
