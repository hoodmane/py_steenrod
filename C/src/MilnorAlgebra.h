//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H
#include <stdbool.h>

#include "combinatorics.h"
#include "Algebra.h"

typedef struct {
    bool generic;
    bool restricted;
    bool truncated;
    uint q_part;
    uint p_part_length;
    uint *p_part;
} Profile;

typedef struct {
    uint q_degree;
    uint q_part;
    uint p_degree;
    uint p_length;
    uint *p_part;
} MilnorBasisElement;

typedef struct {
    uint length;
    MilnorBasisElement *list;
} MilnorBasisElement_list;

typedef struct {
    Algebra algebra;
    bool generic;
    Profile profile;
    int (*sort_order)(const void *a, const void *b); 
} MilnorAlgebra;

Profile *Profile_construct(bool generic, uint q_part_length, uint * q_part, uint p_part_length, uint *p_part, bool truncated);
void Profile_free(Profile *profile);

// Implemented in milnor_datatypes.c
// These methods write a string to a buffer and return the length of the string written.
uint MilnorAlgebra_element_toString(Algebra *algebra, char *buffer, int degree, Vector *m);
void MilnorAlgebra_element_print(char *fmt_string, MilnorAlgebra *algebra, int degree, Vector *m);
uint MilnorMatrix_toString(char *buffer, uint M[MAX_XI_TAU][MAX_XI_TAU], uint rows, uint cols);
uint MilnorAlgebra_basisElement_toKey(char *buffer, MilnorBasisElement *b);
uint MilnorAlgebra_basisElementIndex_toString(Algebra *this, char *buffer, int degree, uint index);
uint MilnorAlgebra_basisElement_toString(char *buffer, MilnorAlgebra *algebra, MilnorBasisElement *b);
char *MilnorAlgebra_basisElement_constructString(MilnorAlgebra *algebra, MilnorBasisElement *b);
void MilnorAlgebra_basisElement_print(char *fmt_string, MilnorAlgebra *algebra, MilnorBasisElement *b);
MilnorBasisElement *MilnorAlgebra_basisElement_fromString(MilnorAlgebra *algebra, char *elt_string);

// Implemented in milnor.c

MilnorAlgebra *MilnorAlgebra_construct(uint p, bool generic, Profile *profile);
void MilnorAlgebra_free(MilnorAlgebra *);

void MilnorAlgebra_generateBasis(Algebra *this, int max_degree);
void MilnorAlgebra_freeBasis(MilnorAlgebra *algebra);

uint MilnorAlgebra_getDimension(Algebra *this, int degree, int excess);
MilnorBasisElement_list MilnorAlgebra_getBasis(MilnorAlgebra *algebra, int degree);
MilnorBasisElement *MilnorAlgebra_basisElement_fromIndex(MilnorAlgebra *algebra, int degree, uint index);
uint MilnorAlgebra_basisElement_toIndex(MilnorAlgebra *algebra,  MilnorBasisElement *b);

void MilnorAlgebra_multiply(Algebra *this, Vector *result, uint coeff, int r_degree, uint r_index, int s_degree, uint s_index, int excess);

#endif //CSTEENROD_MILNOR_H
