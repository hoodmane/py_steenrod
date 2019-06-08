//
// Created by Hood on 5/8/2019.
//

// This file takes care of most of the prime specific basic math.
// It computes mod p binomial and multinomial coefficients, the degrees of taus and xi's,
// and inverses mod p. The function initializePrime(p) must be called before any work
// can be done (without immediately segfaulting).

#ifndef CSTEENROD_COMBINATORICS_H
#define CSTEENROD_COMBINATORICS_H
#include <stdbool.h>

#define MAX_PRIME 251
#define MAX_PRIME_INDEX 54


typedef unsigned int uint;

// Generated with Mathematica:
//   Boole[PrimeQ[#]] PrimePi[#] - 1 & /@ Range[0, 255]
extern int prime_to_index_map[256];

void initializePrime(uint p);
void freePrimes();


int ModPositive(int n, int p);      
uint MinusOneToTheN(uint p, uint n);
uint integer_power(uint a, uint b); // Compute a^b
int power_mod(int p, int b, int e); // Compute b^e mod p

uint logp(uint p, uint n);
void basepExpansion(uint * buffer, uint p, uint n); // Used for Lucas's theorem

int inverse(uint p, int k);

uint Multinomial(uint p, uint len, uint l[]);
uint Multinomial2(uint len, uint l[]);
uint MultinomialOdd(uint p, uint len, uint l[]);

uint Binomial(uint p, int n, int k);
uint Binomial2(uint n, uint k);
uint BinomialOdd(uint p, int n, int k);

#define MAX_XI_TAU 10

uint* getXiDegrees(uint);
uint* getTauDegrees(uint);

#endif //CSTEENROD_COMBINATORICS_H
