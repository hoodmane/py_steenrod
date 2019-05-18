//
// Created by Hood on 4/16/2019.
//

#ifndef CSTEENROD_VECTOR_H
#define CSTEENROD_VECTOR_H

#include <map>
#include <string>
using namespace std;


typedef struct {
    unsigned long degree;
    unsigned long index;
} MonomialIndex;

std::string to_string(MonomialIndex idx);
std::string MonomialIndexToString(const MonomialIndex& idx);

struct cmpMonomialIndex {
    bool operator()(const MonomialIndex& a, const MonomialIndex& b) const;
};

class FpVector {
public:
    long p;
    map<MonomialIndex, long, cmpMonomialIndex> dictionary;
    FpVector(unsigned long);
    void add(FpVector, long);
    void scale(long);
    void add_basis_vector(MonomialIndex, long);
    void assign(FpVector source);
    string toString();

    static FpVector new_basis_vector(unsigned long, MonomialIndex);
};

#endif //CSTEENROD_VECTOR_H
