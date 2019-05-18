//
// Created by Hood on 5/7/2019.
//

#include "FpVector.h"
#include "combinatorics.h"
#include <map>
#include <iostream>

using namespace std;



string to_string(MonomialIndex idx){
    string result = "(";
    result += to_string(idx.degree) + ", ";
    result += to_string(idx.index) + ")";
    return result;
}

string MonomialIndexToString(const MonomialIndex& idx){
    string result = "(";
    result += to_string(idx.degree) + ", ";
    result += to_string(idx.index) + ")";
    return result;
}


bool cmpMonomialIndex::operator()(const MonomialIndex& a, const MonomialIndex& b) const{
    return a.index < b.index;
}


FpVector::FpVector(unsigned long prime) {
    p = prime;
}

void FpVector::add(FpVector v, long c) {
    for(auto it = v.dictionary.begin(); it != v.dictionary.end(); it++){
        this->dictionary[it->first] += c * it->second;
        this->dictionary[it->first] = this->dictionary[it->first] % p;
        if(this->dictionary[it->first] == 0){
            this->dictionary.erase(it->first);
        }
    }
}

void FpVector::add_basis_vector(MonomialIndex idx, long c) {
    dictionary[idx] += c;
    dictionary[idx] = ModPositive(dictionary[idx], p);
    if(dictionary[idx] == 0){
        dictionary.erase(idx);
    }
}

void FpVector::assign(FpVector source) {
    for(auto it = this->dictionary.begin(); it != this->dictionary.end(); it++){
        this->dictionary.erase(it -> first);
    }
    for(auto it = source.dictionary.begin(); it != source.dictionary.end(); it++){
        this->dictionary[it->first] = it->second;
    }
}

void FpVector::scale(long c) {
    for(auto it = this->dictionary.begin(); it != this->dictionary.end(); it++){
        this->dictionary[it->first] *= c;
        this->dictionary[it->first] = ModPositive(this->dictionary[it->first], p);
    }
}

string FpVector::toString() {
    string result = "{";
    for(auto it = this->dictionary.begin(); it != this->dictionary.end(); it++){
        result += MonomialIndexToString(it->first) + ":" + to_string(it->second) + ", ";
    }
    result += "}";
    return result;
}

FpVector FpVector::new_basis_vector(unsigned long prime, MonomialIndex idx){
    FpVector v = FpVector(prime);
    v.dictionary[idx] = 1;
    return v;
}

//int main(){
//    MonomialIndex a, b, c;
//    a.degree = 1;
//    b.degree = 1;
//    c.degree = 1;
//    a.index = 0;
//    b.index = 1;
//    c.index = 2;
//    FpVector v = FpVector::new_basis_vector(5, a);
//    FpVector w = FpVector(5);
//    w.dictionary[a] = 2;
//    w.dictionary[b] = 1;
//    v.add(w, 1);
//    cout << v.toString() << "\n";
//    v.scale(2);
//    cout << v.toString() << "\n";
//    v.scale(3);
//    cout << v.toString() << "\n";
//}