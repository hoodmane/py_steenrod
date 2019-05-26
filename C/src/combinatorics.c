//
// Created by Hood on 4/29/2019.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "combinatorics.h"
#include "FpVector.h"

// Private functions
void initializeInverseTable(uint p);
void initializeBinomialTable(uint p);
void initializeXiTauDegrees(uint p);

uint directBinomial(uint p, uint n, uint k);
uint Multinomial2(uint len, uint* l);
uint Binomial2(uint n, uint k );
uint MultinomialOdd(uint p, uint len, uint* l);
uint BinomialOdd(uint p, uint n, uint k);


int prime_to_index_map[256] = {
    -1, -1, 0, 1, -1, 2, -1, 3, -1, -1, -1, 4, -1, 5, -1, -1, -1, 6, -1, 
    7, -1, -1, -1, 8, -1, -1, -1, -1, -1, 9, -1, 10, -1, -1, -1, -1, -1, 
    11, -1, -1, -1, 12, -1, 13, -1, -1, -1, 14, -1, -1, -1, -1, -1, 15, 
    -1, -1, -1, -1, -1, 16, -1, 17, -1, -1, -1, -1, -1, 18, -1, -1, -1, 
    19, -1, 20, -1, -1, -1, -1, -1, 21, -1, -1, -1, 22, -1, -1, -1, -1, 
    -1, 23, -1, -1, -1, -1, -1, -1, -1, 24, -1, -1, -1, 25, -1, 26, -1, 
    -1, -1, 27, -1, 28, -1, -1, -1, 29, -1, -1, -1, -1, -1, -1, -1, -1, 
    -1, -1, -1, -1, -1, 30, -1, -1, -1, 31, -1, -1, -1, -1, -1, 32, -1, 
    33, -1, -1, -1, -1, -1, -1, -1, -1, -1, 34, -1, 35, -1, -1, -1, -1, 
    -1, 36, -1, -1, -1, -1, -1, 37, -1, -1, -1, 38, -1, -1, -1, -1, -1, 
    39, -1, -1, -1, -1, -1, 40, -1, 41, -1, -1, -1, -1, -1, -1, -1, -1, 
    -1, 42, -1, 43, -1, -1, -1, 44, -1, 45, -1, -1, -1, -1, -1, -1, -1, 
    -1, -1, -1, -1, 46, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 47, 
    -1, -1, -1, 48, -1, 49, -1, -1, -1, 50, -1, -1, -1, -1, -1, 51, -1, 
    52, -1, -1, -1, -1, -1, -1, -1, -1, -1, 53, -1, -1, -1, -1
};

int ModPositive(int  n, int p){
    return ((n % p) + p) % p;
}

uint MinusOneToTheN(uint p, uint n){
    return (n & 1) ? p-1 : 1;
}

uint integer_power(uint b, uint e){
    uint result = 1;
    while(e > 0){
        if((e&1) == 1){
            result *= b;
        }
        b *= b;
        e >>= 1;
    }
    return result;
}

int power_mod(int p, int b, int e){
    int result = 1;
//      b is b^{2^i} mod p
//      if the current bit of e is odd, mutliply b^{2^i} mod p into r.
    while (e > 0){
        if ((e&1) == 1){
            result = (result*b)%p;
        }
        b = (b*b)%p;
        e >>= 1;
    }
    return result;
}

uint *inverse_table[MAX_PRIME_INDEX] = {0};

void initializeInverseTable(uint p){
    uint* table = malloc(p*sizeof(uint));
    for(uint n = 0; n < p; n ++){
        table[n] = power_mod(p, n, p - 2);
    }
    inverse_table[prime_to_index_map[p]] = table;
}


/**
 * Finds the inverse of k mod p.
 * Uses Fermat's little theorem: x^(p-1) = 1 mod p ==> x^(p-2) = x^(-1).
 * @param k an integer
 * @return the inverse of k mod p.
 */
int inverse(uint p, int k){
    return inverse_table[prime_to_index_map[p]][k];
}

uint p_to_the_n_minus_1_over_p_minus_1(uint p, uint n){
    return (integer_power(p, n) - 1) / (p - 1);
}

uint logp(uint p, uint n) {
    uint result = 0;
    while(n > 0){
        n /= p;
        result ++;
    }
    return result;
}

void basepExpansion(uint * result, uint p, uint n){
    uint i = 0;
    for( ; n > 0; n /= p){
        result[i] = n % p;
        i++;
    }
}

uint **binomial_table[MAX_PRIME_INDEX] = {0};

void initializePrime(uint p){
    if(p > MAX_PRIME){
        // Fatal error: max prime allowed is 251
    }
    if(prime_to_index_map[p]==-1){
        // Fatal error: p is not prime.
    }
    if(binomial_table[prime_to_index_map[p]] != NULL){
        return; // Prime already initialized
    }
    initializeBinomialTable(p);
    initializeInverseTable(p);
    initializeXiTauDegrees(p);
    initializeModpLookupTable(p);
    initializeLimbBitIndexLookupTable(p);
}
void freePrimes() {
//    freeBinomialTables();
//    freeInverseTables();
//    freeXiTauDegreess();
}


void initializeBinomialTable(uint p){
    uint** table = (uint**) malloc(p*p*sizeof(uint) + p*sizeof(uint*));
    uint * current_row_ptr = (uint*)(table + p);
    for(uint i = 0; i < p; i++){
        table[i] = current_row_ptr;
        current_row_ptr += p;
    }
    // assert(current_row_ptr == (uint*)(table + p) + p*p);
    for(uint n = 0; n < p; n ++){
        uint entry = 1;
        table[n][0] = entry;
        for(uint k = 1; k <= n; k++) {
            entry *= (n + 1 - k);
            entry /= k;
            table[n][k] = entry % p;
        }
        memset(table[n] + n + 1, 0, (p - n - 1) * sizeof(uint));
    }
    binomial_table[prime_to_index_map[p]] = table;
}

// This is a table lookup, n, k < p.
uint directBinomial(uint p, uint n, uint k){
    return binomial_table[prime_to_index_map[p]][n][k];
}

//Multinomial coefficient of the list l
uint Multinomial2(uint len, uint* l){
    uint bit_or = 0;
    uint sum = 0;
    for(uint i = 0; i < len; i++){
        sum += l[i];
        bit_or |= l[i];
//        if(bit_or < sum){
//            return 0;
//        }
    }
    return (bit_or == sum) ? 1 : 0;
}

//Mod 2 binomial coefficient n choose k
uint Binomial2(uint n, uint k ) {
    if(n < k || k < 0){
        return 0;
    } else {
        if(((n-k) & k) == 0){
            return 1;
        } else {
            return 0;
        }
    }
}

//Mod p multinomial coefficient of l. If p is 2, more efficient to use Multinomial2.
uint MultinomialOdd(uint p, uint len, uint* l){
    uint total = 0;
    for(uint i = 0; i < len; i++){
        total += l[i];
    }
    uint answer = 1;
    uint base_p_expansion_length = logp(p, total);
    uint total_expansion[base_p_expansion_length];
    basepExpansion(total_expansion, p, total);
    uint l_expansions[len][base_p_expansion_length];
    for(uint i=0; i < len; i++){
        memset(l_expansions[i], 0, base_p_expansion_length * sizeof(uint));
        basepExpansion(l_expansions[i], p,  l[i]);
    }
    for(uint index = 0; index < base_p_expansion_length; index++){
        uint multi = 1;
        uint partial_sum = 0;
        for(uint i = 0; i < len; i++){
            partial_sum += l_expansions[i][index];
            if(partial_sum > total_expansion[index]){
                return 0;
            }
            multi *= directBinomial(p, partial_sum, l_expansions[i][index]);
            multi = multi % p;
        }
        answer = (answer * multi) % p;
    }
    return answer;
}

//Mod p binomial coefficient n choose k. If p is 2, more efficient to use Binomial2.
uint BinomialOdd(uint p, uint n, uint k) {
    if( n < k || k < 0 ){
        return 0;
    }
    uint l[2] = { n-k, k };
    return MultinomialOdd(p, 2, l);
}

//Dispatch to Multinomial2 or MultinomialOdd
uint Multinomial(uint p, uint len, uint l[]) {
    if(p == 2){
        return Multinomial2(len, l);
    } else {
        return MultinomialOdd(p, len, l);
    }
}

//Dispatch to Binomial2 or BinomialOdd
uint Binomial(uint p, uint n, uint k){
    if(p == 2){
        return Binomial2(n, k);
    } else {
        return BinomialOdd(p, n, k);
    }
}


uint * xi_degrees[MAX_PRIME_INDEX] = {0};
uint * tau_degrees[MAX_PRIME_INDEX] = {0};

void initializeXiTauDegrees(uint p){
    uint * xi = (uint*)malloc(MAX_XI_TAU * sizeof(uint));
    uint * tau = (uint*)malloc(MAX_XI_TAU * sizeof(uint));
    uint current_xi_degree = 0;
    uint p_to_the_i = 1;
    for(uint i = 0; i < MAX_XI_TAU; i++ ){
        current_xi_degree += p_to_the_i;
        xi[i] = current_xi_degree;
        tau[i] = 2 * p_to_the_i - 1;
        p_to_the_i *= p;
    }
    xi_degrees[prime_to_index_map[p]] = xi;
    tau_degrees[prime_to_index_map[p]] = tau;
}


uint* getTauDegrees(uint p) {
    return tau_degrees[prime_to_index_map[p]];
}

uint* getXiDegrees(uint p) {
    return xi_degrees[prime_to_index_map[p]];
}

//#define C_sdim 5
//#define C_tdim 8
/**
int main(){
    uint p = 7;
    initializePrime(p);
    for(int i=0; i < p; i++){
        printf("%ld -> %ld\n", i, inverse(7, i));
    }
    int src[C_sdim][C_sdim + C_tdim] =
            {{1, 1, 3, 4, 5, 2, 5, 0, 1, 0, 0, 0, 0}, {0, 6, 3, 6, 1, 1, 4, 3, 0, 1, 0, 0, 0}, {3, 0, 3, 3, 0, 2, 6, 3, 0, 0, 1, 0, 0}, {1, 3, 1, 5, 6, 4, 2, 2, 0, 0, 0, 1, 0}, {3, 0, 5, 1, 5, 2, 2, 3, 0, 0, 0, 0, 1}};

    uint rows = C_sdim;
    uint cols = C_sdim + C_tdim;
    int ** M = allocate_matrix(rows, cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            M[i][j] = src[i][j];
        }
    }
    int column_to_pivot_row[cols];
    memset(column_to_pivot_row, 0, cols * sizeof(int));
    row_reduce(M, column_to_pivot_row, p, rows, cols);
    char buffer[200];
    array_to_string(buffer, column_to_pivot_row, cols);
    printf("column_to_pivot_row: %s\n", buffer);
    free(M);
}
**/
