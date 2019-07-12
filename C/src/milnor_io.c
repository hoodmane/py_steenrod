//
// Created by Hood on 5/17/2019.
//

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "combinatorics.h"
#include "FpVector.h"
#include "MilnorAlgebra.h"


// Private functions
uint get_profile_name(char *buffer, Profile P);


uint Profile__getExponent(Profile P, uint p, uint index){
    if(index < P.p_part_length){
        return integer_power(p, P.p_part[index]);
    }
    if(P.truncated){
        return 1;
    }
    return -1;
}

uint Profile_getName(char *buffer, Profile P){
    uint len = 0;
    len += sprintf(buffer + len, "Profile( truncated=%s, ", P.truncated ? "true" : "false");
    if(P.generic){
        len += sprintf(buffer + len, "q_part=%x, ", P.q_part);
    }
    len += sprintf(buffer + len, "p_part=");
    len += array_toString(buffer + len, P.p_part, P.p_part_length);
    len += sprintf(buffer + len, ")");
    return len;
}

Profile *Profile_construct(bool generic, uint q_part_length, uint * q_part, uint p_part_length, uint *p_part, bool truncated){
    Profile * result = malloc(sizeof(Profile) + p_part_length * sizeof(uint));
    result->generic = generic;
    if(q_part_length == 0 && p_part_length == 0 && !truncated){
        result->restricted = false;
        result->truncated = false;
        result->q_part = -1;
        result->p_part_length = 0;
        result->p_part = NULL;  
        return result;      
    }
    result->restricted = true;
    result->truncated = truncated;
    result->generic = generic;
    for(uint i = 0; i<q_part_length; i++){
        result->q_part += (q_part[i] & 1) << i;
    }
    result->p_part_length = p_part_length;
    if(p_part_length > 0){
        result->p_part = (uint*)(result + 1);
        memcpy(result->p_part, p_part, p_part_length*sizeof(uint));
    } else {
        result->p_part = NULL;
    }
    return result;
}

void Profile_free(Profile *profile){
    free(profile);
}

void MilnorAlgebra__generateName(MilnorAlgebra *A){
    // We never initialized it to NULL so this causes trouble.
    // if (A->algebra.name != NULL) {
    //     return;
    // }
    char buffer[200];
    uint len = 0;
    len += sprintf(buffer + len, "MilnorAlgebra(p=%d", A->algebra.p);
    if(A->generic != (A->algebra.p != 2)){
        len += sprintf(buffer + len, ", generic=%s", A->generic ? "true" : "false");
    }
    if(A->profile.restricted){
        len += sprintf(buffer + len, ", ");
        len += Profile_getName(buffer + len, A->profile);
    }
    len += sprintf(buffer + len, ")");
    char *result = malloc((len + 1)* sizeof(char));
    memcpy(result, buffer, len + 1);
    A->algebra.name = result;
}




uint MilnorAlgebra_basisElement_toKey(char *buffer, MilnorBasisElement *b){
    // Copy bytes representing MilnorBasisElement to a string.
//    printf("making hash.\n");
//    printf("    pointer address: %x\n", (int)b->p_part);
    uint len = 0;
    memcpy(buffer + len, &b->q_part, sizeof(uint));
    len += sizeof(uint);
    memcpy(buffer + len, b->p_part, b->p_length * sizeof(uint));
    len += b->p_length * sizeof(uint);

    // Ensure that none of the chars in our key are 0.
    // A zero would terminate our key early.
    for(uint i = 0; i < len; i++){
        buffer[i] = (buffer[i] << 1) | 1;
    }
    // Now add string terminating null character.
    buffer[len] = '\0';
    return len;
}

uint MilnorAlgebra_basisElementIndex_toString(Algebra *this, char *buffer, int degree, uint idx){
    MilnorAlgebra *A = (MilnorAlgebra *)this;
    MilnorBasisElement *b = MilnorAlgebra_basisElement_fromIndex(A, degree, idx);
    return MilnorAlgebra_basisElement_toString(buffer, A, b);
}

uint MilnorAlgebra_basisElement_toString(char *buffer, MilnorAlgebra *A, MilnorBasisElement *b){
    if(b->p_length == 0 && b->q_part == 0){
        buffer[0] = '0';
        buffer[1] = '\0';
        return 1;
    }
    int len = 0;
    char P_or_Sq[4];
    if(A->generic){
        if(b->q_part != 0){
            len += sprintf(buffer + len, "Q(");
            uint idx = 0;
            for(uint q_part = b->q_part; q_part > 0; q_part >>= 1){
                if((q_part & 1) != 0){
                    len += sprintf(buffer + len, "%d,", idx);
                }
                idx ++;
            }
            len += sprintf(buffer + len, ")");
        }
        if(b->p_length != 0 && b->q_part != 0){
            len += sprintf(buffer + len, " * ");
        }
        sprintf(P_or_Sq, "P(");
    } else {
        sprintf(P_or_Sq, "Sq(");
    }
    if(b->p_length != 0){
        len += sprintf(buffer + len, "%s", P_or_Sq);
        for(uint i = 0; i < b->p_length; i++){
            len += sprintf(buffer + len, "%d,", b->p_part[i]);
        }
        len--; // drop trailing comma
        len += sprintf(buffer + len, ")");
    }
    return len;
}

void MilnorAlgebra_basisElement_print(char *fmt_string, MilnorAlgebra *algebra, MilnorBasisElement *b){
    char buffer[2000];
    MilnorAlgebra_basisElement_toString(buffer, algebra, b);
    printf(fmt_string, buffer);
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

// Parses elements like P(0,1) or Q(0,1) * P(0,1)
MilnorBasisElement *MilnorAlgebra_basisElement_fromString(MilnorAlgebra * algebra, char* elt_string){
    char *idx_string, *string, *tofree;
    char * Q_or_P;
    uint p = algebra->algebra.p;
    uint q = algebra->generic ? 2*p-2 : 1;
    tofree = string = strdup(elt_string);

    uint Q_part[20] = {0};
    uint P_part[20] = {0};
    uint Q_part_len = 0;
    uint P_part_len = 0;
    uint q_degree = 0;
    uint p_degree = 0;

    while ((Q_or_P = strsep(&string, "*")) != NULL){
        Q_or_P = trimwhitespace(Q_or_P);
        uint* target = Q_part;
        uint len = 0;
        uint *Q_or_P_len = &Q_part_len;
        if(Q_or_P[0] == 'P'){
            target = P_part;
            Q_or_P_len = &P_part_len;
        }
        if(Q_or_P[0] == 'S'){// Presumably the S is part of "Sq"
            target = P_part;
            Q_or_P_len = &P_part_len;
            Q_or_P++;// "Sq("" is three characters, one longer than "Q(" or "P("
        }
        Q_or_P += 2; // "Q(" and "P(" are two characters
        Q_or_P[strlen(Q_or_P)-1] = '\0';
        char *end;
        while((idx_string = strsep(&Q_or_P, ",")) != NULL){
            target[len] = strtoul(idx_string, &end, 10);
            len ++;
        }
        *Q_or_P_len = len;
    }
    free(tofree);
    tofree = NULL;
    string = NULL;
    int* tau_degrees = getTauDegrees(p);
    int* xi_degrees = getXiDegrees(p);
    uint Q_bit_string = 0;
    for(uint i = 0; i < Q_part_len; i++){
        Q_bit_string += (1 << Q_part[i]);
        q_degree += tau_degrees[Q_part[i]];
    }

    for(uint i = 0; i < P_part_len; i++){
        p_degree += xi_degrees[i] * P_part[i] * q;
    }

    MilnorBasisElement result;
    result.p_part = P_part;
    result.p_length = P_part_len;
    result.p_degree = p_degree;
    result.q_part = Q_bit_string;
    result.q_degree = q_degree;

    uint idx = MilnorAlgebra_basisElement_toIndex(algebra, &result);
    return MilnorAlgebra_basisElement_fromIndex(algebra, result.p_degree + result.q_degree, idx);
}

uint MilnorAlgebra_element_toString(Algebra *this, char *buffer, int degree, Vector * m){
    MilnorAlgebra *algebra = (MilnorAlgebra*)this;
    uint len = 0;
    for(
        VectorIterator it = Vector_getIterator(m); 
        it.has_more; 
        it = Vector_stepIterator(it)
    ){
        if(it.value == 0){
            continue;
        }
        if(it.value != 1) {
            len += sprintf(buffer + len, "%d * ", it.value);
        }
        MilnorBasisElement *b = MilnorAlgebra_basisElement_fromIndex(algebra, degree, it.index);
        len += MilnorAlgebra_basisElement_toString(buffer + len, algebra, b);
        len += sprintf(buffer + len, " + ");
    }
    if(len == 0){
        len += sprintf(buffer + len, "0");
    } else {
        // Remove trailing " + "
        len -= 3;
        buffer[len] = '\0';
    }
    return len;
}

void MilnorElement_print(char *fmt_string, MilnorAlgebra *algebra, uint degree, Vector *m){
    char buffer[2000];
    uint len = MilnorAlgebra_element_toString((Algebra *)algebra, buffer, degree, m);
    assert(len < 2000);
    printf(fmt_string, buffer);
}

uint MilnorMatrix_toString(char *buffer, uint M[MAX_XI_TAU][MAX_XI_TAU], uint rows, uint cols){
    uint len = 0;
    len += sprintf(buffer + len, "[\n");
    for(uint row = 0; row < rows; row++) {
        len += sprintf(buffer + len, "  ");
        len += array_toString(buffer + len, M[row], cols);
        len += sprintf(buffer + len, ",");
    }
    len += sprintf(buffer + len, "]\n");
    return len;
}