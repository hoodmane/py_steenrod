//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_COMBINATORICS_H
#define CSTEENROD_COMBINATORICS_H
#include <stdbool.h>

#define MAX_PRIME 251
#define MAX_PRIME_INDEX 54


typedef unsigned int uint;

// Generated with Mathematica:
//   Boole[PrimeQ[#]] PrimePi[#] - 1 & /@ Range[0, 255]
extern int prime_to_index_map[256];

int ModPositive(int n, int p);
uint MinusOneToTheN(uint p, uint n);
uint integer_power(uint a, uint b);
int power_mod(int p, int b, int e);

uint p_to_the_n_minus_1_over_p_minus_1(uint p, uint n);
uint logp(uint p, uint n);
void basepExpansion(uint * buffer, uint p, uint n);

void initializePrime(uint p);
void freePrimes();

int inverse(uint p, int k);

uint Multinomial(uint p, uint len, uint l[]);
uint Binomial(uint p, uint n, uint k);

#define MAX_XI_TAU 10

uint* getXiDegrees(uint);
uint* getTauDegrees(uint);

uint getBitlength(uint);
uint getBitMask(uint p);
uint getEntriesPer64Bits(uint);

uint modPLookup(uint p, uint n);

#endif //CSTEENROD_COMBINATORICS_H
