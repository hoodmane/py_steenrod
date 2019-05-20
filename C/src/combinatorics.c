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
    initializeInverseTable(p);
}

long ModPositive(long  n, long p){
    return ((n % p) + p) % p;
}

long MinusOneToTheN(long n){
    return -(n % 2 * 2 - 1);
}

unsigned long integer_power(unsigned long b, unsigned long e){
    unsigned long result = 1;
    while(e > 0){
        if((e&1) == 1){
            result *= b;
        }
        b *= b;
        e >>= 1;
    }
    return result;
}

long power_mod(long p, long b, long e){
    long result = 1;
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

khash_t(prime_to_list) * inverse_table;

void initializeInverseTable(unsigned long p){
    if(inverse_table == NULL){
        inverse_table = kh_init(prime_to_list);
    }
    khint_t k;
    int absent;
    k = kh_put(prime_to_list, inverse_table, p, &absent);  // insert a key to the hash table
    if(!absent){
        return;
    }
    unsigned long* table_p = malloc(p*sizeof(unsigned long*));
    for(unsigned long n = 0; n < p; n ++){
        table_p[n] = power_mod(p, n, p - 2);
    }
    kh_val(inverse_table, k) = table_p;
}


/**
 * Finds the inverse of k mod p.
 * Uses Fermat's little theorem: x^(p-1) = 1 mod p ==> x^(p-2) = x^(-1).
 * @param k an integer
 * @return the inverse of k mod p.
 */
long inverse(unsigned long p, long k){
    khint_t key = kh_get(prime_to_list, inverse_table, p);
    return kh_val(inverse_table, key)[k];
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
        memset(table_p[n] + n + 1, 0, (p - n - 1) * sizeof(unsigned long));
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

long** allocate_matrix(unsigned long rows, unsigned long cols)  {
    long** M = calloc(1, rows * sizeof(long*) + rows * cols * sizeof(long));
    for(long row = 0; row < rows; row++){
        M[row] = (long*)(M + rows) + row * cols;
    }
    return M;
}

void print_matrix(long **M, unsigned long rows, unsigned long columns){
    printf("    [\n");
    for(long i = 0; i < rows; i++){
        printf("        [");
        for(long j = 0; j < columns; j++){
            printf("%ld, ", M[i][j]);
        }
        printf("]\n");
    }
    printf("    ]\n");
}

void row_reduce(row_reduce_state * state){
    unsigned long p = state->p;
    unsigned long source_dim = state->source_dim;
    unsigned long target_dim = state->target_dim;
    long ** matrix = state->matrix;
    for(state->pivot++ ; state->pivot < target_dim; state->pivot++){
        long pivot_row;
        for(pivot_row = state->pivot; pivot_row < source_dim; pivot_row ++){
            if(matrix[pivot_row][state->pivot] != 0){
                break;
            }
        }
        if(pivot_row == source_dim){
            state->found_cokernel = true;
            return;
        }
        print_matrix(matrix, state->source_dim, state->source_dim + state->target_dim);
        long * temp = matrix[state->pivot];
        matrix[state->pivot] = matrix[pivot_row];
        matrix[pivot_row] = temp;
        printf("row(%ld) <==> row(%ld)\n", state->pivot, pivot_row);
        print_matrix(matrix, state->source_dim, state->source_dim + state->target_dim);

        long c = matrix[state->pivot][state->pivot];
        long c_inv = inverse(state->p, c);
        for(long column = state -> pivot; column < source_dim + target_dim; column ++){
            matrix[state->pivot][column] = (matrix[state->pivot][column] * c_inv) % p;
        }
        printf("row(%ld) *= %ld\n", state->pivot, c_inv);
        print_matrix(matrix, state->source_dim, state->source_dim + state->target_dim);
        for(long row = 0; row < state->pivot; row++){
            long row_op_coeff = (-matrix[row][state->pivot] + p) % p;
            if(row_op_coeff == 0){
                continue;
            }
            // Do row operation
            for(long column = state -> pivot; column < source_dim + target_dim; column++){
                matrix[row][column] = (matrix[row][column] + row_op_coeff * matrix[state->pivot][column]) % p;
            }
            printf("row(%ld) += %ld * row(%ld)\n", row, row_op_coeff, state->pivot);
            print_matrix(matrix, state->source_dim, state->source_dim + state->target_dim);
        }
        // Between pivot and pivot_row, we already checked that the pivot column is 0, so skip ahead a bit.
        for(long row = pivot_row + 1; row < source_dim; row++){
            long row_op_coeff = (-matrix[row][state->pivot] + p) % p;
            if(row_op_coeff == 0){
                continue;
            }
            // Do row operation
            for(long column = state -> pivot; column < source_dim + target_dim; column++){
                matrix[row][column] = (matrix[row][column] + row_op_coeff * matrix[state->pivot][column]) % p;
            }
            printf("row(%ld) += %ld * row(%ld)\n", row, row_op_coeff, state->pivot);
            print_matrix(matrix, state->source_dim, state->source_dim + state->target_dim);
        }
    }
    state->found_cokernel = false;
    return;
}



#define C_sdim 5
#define C_tdim 8
/**
int main(){
    unsigned long p = 7;
    initializePrime(p);
    for(long i=0; i < p; i++){
        printf("%ld -> %ld\n", i, inverse(7, i));
    }
    long src[C_sdim][C_sdim + C_tdim] =
            {{1, 1, 3, 4, 5, 2, 5, 0, 1, 0, 0, 0, 0}, {0, 6, 3, 6, 1, 1, 4, 3, 0, 1, 0, 0, 0}, {3, 0, 3, 3, 0, 2, 6, 3, 0, 0, 1, 0, 0}, {1, 3, 1, 5, 6, 4, 2, 2, 0, 0, 0, 1, 0}, {3, 0, 5, 1, 5, 2, 2, 3, 0, 0, 0, 0, 1}};
    row_reduce_state state;
    state.p = p;
    state.row_capacity = 100;
    state.column_capacity = 100;
    state.pivot = -1;
    state.source_dim = C_sdim;
    state.target_dim = C_tdim;
    long ** M = allocate_matrix(100, 100);
    state.matrix = M;
    for(long i = 0; i < state.source_dim; i++){
        for(long j = 0; j < state.target_dim + state.source_dim; j++){
            M[i][j] = src[i][j];
        }
    }
    row_reduce(&state);
    printf("[\n");
    for(long i = 0; i < state.source_dim; i++){
        printf("    [");
        for(long j = 0; j < state.target_dim + state.source_dim; j++){
            printf("%ld, ", M[i][j]);
        }
        printf("]\n");
    }
    printf("]\n");
    printf("found_cokernel? %s\n", state.found_cokernel ? "true" : "false");
    printf("pivot: %ld\n", state.pivot);
    free(state.matrix);
}
/**/