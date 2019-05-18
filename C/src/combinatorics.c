//
// Created by Hood on 4/29/2019.
//

#include "combinatorics.h"
#include <stdio.h>

#include "khash.h"
KHASH_MAP_INIT_INT(prime_to_table, unsigned long**)
KHASH_MAP_INIT_INT(prime_to_list, unsigned long*)

void initializePrime(unsigned long p){
    directBinomialInitializeTable(p);
    initializeXiTauDegrees(p);
}

long ModPositive(long  n, long p){
    return ((n % p) + p) % p;
}

long MinusOneToTheN(long n){
    return -(n % 2 * 2 - 1);
}

unsigned long integer_power(unsigned long a, unsigned long b){
    unsigned long res = 1;
    while(b > 0 ){
        if(b & 1 != 0){
            res *= a;
        }
        a *= a;
        b >>= 1;
    }
    return res;
}

unsigned long p_to_the_n_minus_1_over_p_minus_1(unsigned long p, unsigned long n){
    return (integer_power(p, n) - 1) / (p - 1);
}

unsigned long logp(unsigned long p, unsigned long n) {
    unsigned long result = 0;
    while(n > 0){
        n /= p;
        result ++;
    }
    return result;
}

unsigned long* basepExpansion(unsigned long p, unsigned long n, unsigned long padlength){
    unsigned long* result = (unsigned long*) calloc(padlength, sizeof(unsigned long));
    unsigned long i = 0;
    for( ; n > 0; n /= p){
        result[i] = n % p;
        i++;
    }
    return result;
}

khash_t(prime_to_table) * binomial_table = NULL;

void directBinomialInitializeTable(unsigned long p){
    if(binomial_table == NULL){
        binomial_table = kh_init(prime_to_table);
    }
    khint_t k;
    int absent;
    k = kh_put(prime_to_table, binomial_table, p, &absent);  // insert a key to the hash table
    if(!absent){
        return;
    }
    unsigned long** table_p = (unsigned long**) malloc(p*p*sizeof(unsigned long) + p*sizeof(unsigned long*));
    for(unsigned long i = 0; i < p; i++){
        table_p[i] = (unsigned long*) (table_p + p) + p*i;
    }
    for(unsigned long n = 0; n < p; n ++){
        unsigned long entry = 1;
        table_p[n][0] = entry;
        for(unsigned long k = 1; k <= n; k++) {
            entry *= (n + 1 - k);
            entry /= k;
            table_p[n][k] = entry % p;
        }
    }
    kh_val(binomial_table, k) = table_p;
}

unsigned long directBinomial(unsigned long p, unsigned long n, unsigned long k){
    khint_t key = kh_get(prime_to_table, binomial_table, p);
    return kh_val(binomial_table, key)[n][k];
}

//Multinomial coefficient of the list l
unsigned long Multinomial2(unsigned long len, unsigned long* l){
    unsigned long bit_or = 0;
    unsigned long sum = 0;
    for(unsigned long i = 0; i < len; i++){
        sum += l[i];
        bit_or |= l[i];
//        if(bit_or < sum){
//            return 0;
//        }
    }
    return (bit_or == sum) ? 1 : 0;
}

//Mod 2 binomial coefficient n choose k
unsigned long Binomial2(unsigned long n, unsigned long k ) {
    if(n < k || k < 0){
        return 0;
    } else {
        if((n-k) & k == 0){
            return 1;
        } else {
            return 0;
        }
    }
}

//Mod p multinomial coefficient of l. If p is 2, more efficient to use Multinomial2.
unsigned long MultinomialOdd(unsigned long p, unsigned long len, unsigned long* l){
    unsigned long total = 0;
    for(unsigned long i = 0; i < len; i++){
        total += l[i];
    }
    unsigned long answer = 1;
    unsigned long base_p_expansion_length = logp(p, total);
    unsigned long* total_expansion = basepExpansion(p, total, base_p_expansion_length);
    unsigned long *l_expansions[len];
    for(unsigned long i=0; i < len; i++){
        l_expansions[i] = basepExpansion(p,  l[i], base_p_expansion_length );
    }
    for(unsigned long index = 0; index < base_p_expansion_length; index++){
        unsigned long multi = 1;
        unsigned long partial_sum = 0;
        for(unsigned long i = 0; i < len; i++){
            partial_sum += l_expansions[i][index];
            if(partial_sum > total_expansion[index]){
                answer = 0;
                goto BREAK;
            }
            multi *= directBinomial(p, partial_sum, l_expansions[i][index]);
            multi = multi % p;
        }
        answer = (answer * multi) % p;
    }
    BREAK:
    for(unsigned long i = 0; i < len; i++){
        free(l_expansions[i]);
    }
    return answer;
}

//Mod p binomial coefficient n choose k. If p is 2, more efficient to use Binomial2.
unsigned long BinomialOdd(unsigned long p, unsigned long n, unsigned long k) {
    if( n < k || k < 0 ){
        return 0;
    }
    unsigned long l[2] = { n-k, k };
    return MultinomialOdd(p, 2, l);
}

//Dispatch to Multinomial2 or MultinomialOdd
unsigned long Multinomial(unsigned long p, unsigned long len, unsigned long l[]) {
    if(p == 2){
        return Multinomial2(len, l);
    } else {
        return MultinomialOdd(p, len, l);
    }
}

//Dispatch to Binomial2 or BinomialOdd
unsigned long Binomial(unsigned long p, unsigned long n, unsigned long k){
    if(p == 2){
        return Binomial2(n, k);
    } else {
        return BinomialOdd(p, n, k);
    }
}


khash_t(prime_to_list) * xi_degrees, * tau_degrees;

void initializeXiTauDegrees(unsigned long p){
    if(xi_degrees == NULL){
        xi_degrees = kh_init(prime_to_list);
        tau_degrees = kh_init(prime_to_list);
    }
    khint_t xi_key, tau_key;
    int absent;
    xi_key = kh_put(prime_to_list, xi_degrees, p, &absent);  // insert a key to the hash table
    if(!absent){
        return;
    }
    tau_key = kh_put(prime_to_list, tau_degrees, p, &absent);  // insert a key to the hash table

    kh_val(xi_degrees, xi_key) = (unsigned long*)malloc(MAX_XI_TAU * sizeof(unsigned long));
    kh_val(tau_degrees, tau_key) = (unsigned long*)malloc(MAX_XI_TAU * sizeof(unsigned long));

    unsigned long current_xi_degree = 0;
    unsigned long p_to_the_i = 1;
    for(unsigned long i = 0; i < MAX_XI_TAU; i++ ){
        current_xi_degree += p_to_the_i;
        kh_val(xi_degrees, xi_key)[i] = current_xi_degree;
        kh_val(tau_degrees, tau_key)[i] = 2 * p_to_the_i - 1;
        p_to_the_i *= p;
    }
}


unsigned long* getTauDegrees(unsigned long p) {
    khint_t key;
    key = kh_get(prime_to_list, tau_degrees, p);
//    if(key == kh_end(tau_degrees)){
//        initializeXiTauDegrees(p);
//        key = kh_get(prime_to_list, tau_degrees, p);
//    }
    return kh_value(tau_degrees, key);
}

unsigned long* getXiDegrees(unsigned long p) {
    khint_t key;
    key = kh_get(prime_to_list, xi_degrees, p);
    if(key == kh_end(xi_degrees)){
        initializeXiTauDegrees(p);
        key = kh_get(prime_to_list, xi_degrees, p);
    }
    return kh_value(xi_degrees, key);
}


//int main(){
//    directBinomialInitializeTable(7);
//    for(int n = 0; n < 7; n++){
//        for(int k = 0; k < 7; k++){
//            cout << directBinomial(7, n, k) << " ";
//        }
//        cout << "\n";
//    }
//    int l[5] = {1, 2, 4, 8, 15};
//    cout << Multinomial2(5, l) << "\n";
//
//    map<int, int> test;
//    test[0] += 5;
//    cout << test[0] << "\n";
//
//}