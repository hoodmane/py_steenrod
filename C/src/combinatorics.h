//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_COMBINATORICS_H
#define CSTEENROD_COMBINATORICS_H

void initializePrime(unsigned long p);

long ModPositive(long n, long p);
long MinusOneToTheN(long n);
unsigned long integer_power(unsigned long a, unsigned long b);
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

#endif //CSTEENROD_COMBINATORICS_H
