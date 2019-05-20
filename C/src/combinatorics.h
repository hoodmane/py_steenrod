//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_COMBINATORICS_H
#define CSTEENROD_COMBINATORICS_H
#include <stdbool.h>

void initializePrime(unsigned long p);

long ModPositive(long n, long p);
long MinusOneToTheN(long n);
unsigned long integer_power(unsigned long a, unsigned long b);

long power_mod(long p, long b, long e);
void initializeInverseTable(unsigned long p);
long inverse(unsigned long p, long k);

unsigned long p_to_the_n_minus_1_over_p_minus_1(unsigned long p, unsigned long n);
unsigned long logp(unsigned long p, unsigned long n);
//
unsigned long* basepExpansion(unsigned long p, unsigned long n, unsigned long padlength);

void directBinomialInitializeTable(unsigned long p);
unsigned long directBinomial(unsigned long p, unsigned long n, unsigned long k);
unsigned long Multinomial2(unsigned long len, unsigned long* l);
unsigned long Binomial2(unsigned long n, unsigned long k );
unsigned long MultinomialOdd(unsigned long p, unsigned long len, unsigned long* l);
unsigned long BinomialOdd(unsigned long p, unsigned long n, unsigned long k);
unsigned long Multinomial(unsigned long p, unsigned long len, unsigned long l[]);
unsigned long Binomial(unsigned long p, unsigned long n, unsigned long k);

#define MAX_XI_TAU 10
void initializeXiTauDegrees(unsigned long);
unsigned long* getXiDegrees(unsigned long);
unsigned long* getTauDegrees(unsigned long);


long** allocate_matrix(unsigned long rows, unsigned long cols);

typedef struct {
    unsigned long p;
    long source_dim;
    long target_dim;
    long pivot;
    long row_capacity;
    long column_capacity;
    long ** matrix;
    bool found_cokernel;
} row_reduce_state;

void row_reduce(row_reduce_state * state);

#endif //CSTEENROD_COMBINATORICS_H
