//
// Created by Hood on 5/17/2019.
//

#include <stdio.h>
#include <ctype.h>

#include "milnor_datatypes.h"
#include "milnor.h"
#include "combinatorics.h"

unsigned long getProfileExponent(Profile P, unsigned long p, unsigned long index){
    if(index < P.p_part_length){
        return integer_power(p, P.p_part[index]);
    }
    if(P.truncated){
        return 1;
    }
    return -1;
}

bool checkProfile(unsigned long profile[], P_part xi_mono){
    for(int i = 0; i < xi_mono.length; i++){
        if(xi_mono.p_part[i] >= profile[i]){
            return false;
        }
    }
    return true;
}

void generate_profile_name(Profile P){
    if (P.name != NULL) {
        return;
    }
    char buffer[200];
    unsigned long len = 0;
    len += sprintf(buffer + len, "Profile( truncated=%s, ", P.truncated ? "true" : "false");
    len += sprintf(buffer + len, "q_part=%lx, ", P.q_part);
    len += sprintf(buffer + len, "p_part=");
    len += array_to_string(buffer, P.p_part, P.p_part_length);
    len += sprintf(buffer + len, ")");
    P.name = malloc((len + 1)* sizeof(char));
    memcpy(P.name, buffer, len);
}

void initializeMilnorAlgebraFields(MilnorAlgebra * algebra){
    algebra->profile.truncated = false;
    algebra->profile.p_part_length = 0;
    algebra->profile.q_part = -1;

    algebra -> P_table = NULL;
    algebra -> P_table_by_P_length = NULL;
    algebra->P_table_max_degree = 0;

    algebra->Q_table = NULL;
    algebra->Q_table_max_tau = 0;

    algebra->basis_table = NULL;
    algebra->basis_max_degree = -1;
    algebra->basis_name_to_index_map = NULL;
}

void milnor_algebra_generate_name(MilnorAlgebra *A){
    if (A->name != NULL) {
        return;
    }
    char buffer[200];
    long len = 0;
    len += sprintf(buffer + len, "MilnorAlgebra(p=%ld, generic=%s", A->p, A->generic ? "true" : "false" );
    if(A->profile.restricted){
        generate_profile_name(A->profile);
        len += sprintf(buffer + len, "%s", A->profile.name);
    }
    len += sprintf(buffer + len, ")");
    char * result = malloc((len + 1)* sizeof(char));
    memcpy(result, buffer, len + 1);
    A->name = result;
}


int array_to_string(string buffer, unsigned long* A, unsigned long length){
    buffer[0] = '[';
    buffer[1] = '\0';
    long len = 1;
    for (int i = 0; i < length; i++) {
        len += sprintf(buffer + len, "%ld, ", A[i]);
    }
    len += sprintf(buffer + len, "]");
    return len;
}


int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b){
    if(b->p_length == 0 && b->q_part == 0){
        buffer[0] = '0';
        buffer[1] = '\0';
        return 1;
    }
    long len = 0;
    if(b->q_part != 0){
        len += sprintf(buffer + len, "Q(");
        unsigned long idx = 0;
        for(unsigned long q_part = b->q_part; q_part > 0; q_part >>= 1){
            if(q_part & 1 != 0){
                len += sprintf(buffer + len, "%ld,", idx);
            }
            idx ++;
        }
        len += sprintf(buffer + len, ")");
    }
    if(b->p_length != 0 && b->q_part != 0){
        len += sprintf(buffer + len, " ");
    }
    if(b->p_length != 0){
        len += sprintf(buffer + len, "P(");
        for(unsigned long i = 0; i < b->p_length; i++){
            len += sprintf(buffer + len, "%ld,", b->p_part[i]);
        }
        len += sprintf(buffer + len, ")");
    }
    return len;
}

// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char *trimwhitespace(char *str){
    char *end;

    // Trim leading space
    for(;isspace((unsigned char)*str); str++){ };

    if(*str == 0) {  // All spaces?
        return str;
    }
    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    end[1] = '\0';
    return str;
}

MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string){
    char *idx_string, *string, *tofree;
    char * Q_or_P;
    unsigned long p = algebra->p;
    unsigned long q = algebra->generic ? 2*p-2 : 1;
    tofree = string = strdup(elt_string);

    unsigned long Q_part[20] = {0};
    unsigned long P_part[20] = {0};
    unsigned long Q_part_len = 0;
    unsigned long P_part_len = 0;
    unsigned long q_degree = 0;
    unsigned long p_degree = 0;

    while ((Q_or_P = strsep(&string, "*")) != NULL){
        Q_or_P = trimwhitespace(Q_or_P);
        unsigned long int* target = Q_part;
        unsigned long int len = 0;
        unsigned long int *Q_or_P_len = &Q_part_len;
        if(Q_or_P[0] == 'P'){
            target = P_part;
            Q_or_P_len = &P_part_len;
        }
        Q_or_P += 2;
        Q_or_P[strlen(Q_or_P)-1] = '\0';
        char * end;
        while((idx_string = strsep(&Q_or_P, ",")) != NULL){
            target[len] = strtol(idx_string, &end, 10);
            len ++;
        }
        *Q_or_P_len = len;
    }
    free(tofree);
    tofree = NULL;
    string = NULL;
    unsigned long* tau_degrees = getTauDegrees(p);
    unsigned long* xi_degrees = getXiDegrees(p);
    unsigned long Q_bit_string = 0;
    for(unsigned long i = 0; i < Q_part_len; i++){
        Q_bit_string += (1 << Q_part[i]);
        q_degree += tau_degrees[Q_part[i]];
    }

    unsigned long * output_P_part = (unsigned long*) malloc(P_part_len * sizeof(unsigned long));
    memcpy(output_P_part, P_part, P_part_len * sizeof(unsigned long));
    for(unsigned long i = 0; i < P_part_len; i++){
        p_degree += xi_degrees[i] * P_part[i] * q;
    }

    MilnorBasisElement result;
    result.p_part = output_P_part;
    result.p_length = P_part_len;
    result.p_degree = p_degree;
    result.q_part = Q_bit_string;
    result.q_degree = q_degree;
    return result;
}

int milnor_element_to_string(string buffer, MilnorAlgebra * algebra, Vector * m){
    unsigned long len = 0;
    for(unsigned long i = 0; i < m->dimension; i ++){
        if(m->vector[i] == 0){
            continue;
        }
        if(m->vector[i] != 1) {
            len += sprintf(buffer + len, "%ld * ", m->vector[i]);
        }
        MilnorBasisElement b = GetMilnorBasisElementFromIndex(algebra, m->degree, i);
        len += milnor_basis_element_to_string(buffer + len, &b);
        len += sprintf(buffer + len, " + ");
    }
    if(len == 0){
        len += sprintf(buffer + len, "0");
    }
    return len;
}

int milnor_matrix_to_string(string buffer, unsigned long** M, unsigned long rows, unsigned long cols){
    unsigned long len = 0;
    len += sprintf(buffer + len, "[\n");
    for(int row = 0; row < rows; row++) {
        len += sprintf(buffer + len, "  ");
        len += array_to_string(buffer + len, M[row], cols);
        len += sprintf(buffer + len, ",");
    }
    len += sprintf(buffer + len, "]\n");
    return len;
}