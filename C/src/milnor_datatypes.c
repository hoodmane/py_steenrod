//
// Created by Hood on 5/17/2019.
//

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "milnor.h"
#include "FpVector.h"
#include "combinatorics.h"

// Private functions
uint get_profile_name(char *buffer, Profile P);


uint Profile_getExponent(Profile P, uint p, uint index){
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
    len += sprintf(buffer + len, "q_part=%x, ", P.q_part);
    len += sprintf(buffer + len, "p_part=");
    len += array_toString(buffer, P.p_part, P.p_part_length);
    len += sprintf(buffer + len, ")");
    return len;
}

void MilnorAlgebra_generateName(MilnorAlgebra *A){
    if (A->name != NULL) {
        return;
    }
    char buffer[200];
    int len = 0;
    len += sprintf(buffer + len, "MilnorAlgebra(p=%d, generic=%s", A->p, A->generic ? "true" : "false" );
    if(A->profile.restricted){
        len += Profile_getName(buffer + len, A->profile);
    }
    len += sprintf(buffer + len, ")");
    char * result = malloc((len + 1)* sizeof(char));
    memcpy(result, buffer, len + 1);
    A->name = result;
}





int MilnorBasisElement_toKey(string buffer, MilnorBasisElement *b){
    // Copy bytes representing MilnorBasisElement to a string.
//    printf("making hash.\n");
//    printf("    pointer address: %x\n", (int)b->p_part);
    memcpy(buffer, &b->q_part, sizeof(uint));
    memcpy(buffer + sizeof(uint), b->p_part, b->p_length * sizeof(uint));
    int len = ((b->p_length + 1) * sizeof(uint)) / sizeof(char);

    // Now we need to make sure that none of the entries are 0.
    // That would end the string early.
    for(unsigned int i = 0; i < len; i++){
        buffer[i]++;
    }
    // Now add string terminating null character.
    buffer[len] = '\0';
    return len;
}

int MilnorBasisElement_toString(string buffer, MilnorBasisElement *b){
    if(b->p_length == 0 && b->q_part == 0){
        buffer[0] = '0';
        buffer[1] = '\0';
        return 1;
    }
    int len = 0;
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
        len += sprintf(buffer + len, " ");
    }
    if(b->p_length != 0){
        len += sprintf(buffer + len, "P(");
        for(uint i = 0; i < b->p_length; i++){
            len += sprintf(buffer + len, "%d,", b->p_part[i]);
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

MilnorBasisElement MilnorBasisElement_fromString(MilnorAlgebra * algebra, char* elt_string){
    char *idx_string, *string, *tofree;
    char * Q_or_P;
    uint p = algebra->p;
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
    uint* tau_degrees = getTauDegrees(p);
    uint* xi_degrees = getXiDegrees(p);
    uint Q_bit_string = 0;
    for(uint i = 0; i < Q_part_len; i++){
        Q_bit_string += (1 << Q_part[i]);
        q_degree += tau_degrees[Q_part[i]];
    }

    uint * output_P_part = (uint*) malloc(P_part_len * sizeof(uint));
    memcpy(output_P_part, P_part, P_part_len * sizeof(uint));
    for(uint i = 0; i < P_part_len; i++){
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

int MilnorElement_toString(string buffer, MilnorAlgebra * algebra, uint degree, Vector * m){
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
        MilnorBasisElement b = MilnorBasisElement_fromIndex(algebra, degree, it.index);
        len += MilnorBasisElement_toString(buffer + len, &b);
        len += sprintf(buffer + len, " + ");
    }
    if(len == 0){
        len += sprintf(buffer + len, "0");
    }
    return len;
}

int milnor_matrix_to_string(string buffer, uint M[MAX_XI_TAU][MAX_XI_TAU], uint rows, uint cols){
    uint len = 0;
    len += sprintf(buffer + len, "[\n");
    for(int row = 0; row < rows; row++) {
        len += sprintf(buffer + len, "  ");
        len += array_toString(buffer + len, M[row], cols);
        len += sprintf(buffer + len, ",");
    }
    len += sprintf(buffer + len, "]\n");
    return len;
}