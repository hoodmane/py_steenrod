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

// This computes all the tables relevant for the functions in combinatorics.h and the ones in FpVector.h.
// This must be called before doing any work at the prime p or else we segfault.
void initializePrime(uint p);
void freePrimes();

// Compute n mod p in the range from 0 to p-1. Pretty much unnecessary?
int ModPositive(int n, int p);      
// Compute (-1)^n mod p. It either returns 1 if n is even or p-1 if n is odd.
uint MinusOneToTheN(uint p, uint n);
// Compute a^b
uint integer_power(uint a, uint b); 
// Compute b^e mod p
int power_mod(int p, int b, int e); 

// Compute ceil(logp(n)) -- so the length of the base p expansion of n
uint logp(uint p, uint n);

// Compute the base p expansion of n, write result into buffer. Used for Lucas's theorem.
void basepExpansion(uint * buffer, uint p, uint n); 

// Compute k^{-1} mod p.
int inverse(uint p, int k);

// Compute the mod p multinomial coefficient of a list l of length len.
uint Multinomial(uint p, uint len, uint l[]);
uint Multinomial2(uint len, uint l[]);
uint MultinomialOdd(uint p, uint len, uint l[]);

// Compute the mod p binomial coefficient n choose k.
uint Binomial(uint p, int n, int k);
uint Binomial2(uint n, uint k);
uint BinomialOdd(uint p, int n, int k);

#define MAX_XI_TAU 10

// Get the degrees of xi_i and tau_i. 
// We only work with xi_i for i < 10 -- xi_10 is at least of degree 2^10 = 1024, 
// so probably we'll never get out that far.
int* getXiDegrees(uint);
int* getTauDegrees(uint);

#endif //CSTEENROD_COMBINATORICS_H
