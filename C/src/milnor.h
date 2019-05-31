//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H
#include <stdbool.h>

#include "combinatorics.h"
#include "algebra.h"

typedef char *string;

typedef struct {
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
    uint p;
    bool generic;
    Profile profile;
    string name;
    uint max_degree;
} MilnorAlgebra;


// Implemented in milnor_datatypes.c
// These methods write a string to a buffer and return the length of the string written.
int MilnorElement_toString(string buffer, MilnorAlgebra *algebra, uint degree, Vector *m);
int MilnorMatrix_toString(string buffer, uint M[MAX_XI_TAU][MAX_XI_TAU], uint rows, uint cols);
int MilnorBasisElement_toKey(string buffer, MilnorBasisElement *b);
int MilnorBasisElement_toString(string buffer, MilnorBasisElement *b);
MilnorBasisElement MilnorBasisElement_fromString(MilnorAlgebra *algebra, char *elt_string);

MilnorAlgebra *MilnorAlgebra_construct(uint p, bool generic, Profile *profile);
void MilnorAlgebra_free(MilnorAlgebra *);

// Implemented in milnor.c
bool MilnorAlgebra_generateBasis(Algebra *algebra, uint max_degree);
void MilnorAlgebra_freeBasis(MilnorAlgebra *algebra);

uint MilnorAlgebra_getDimension(Algebra *algebra, uint degree);
MilnorBasisElement_list MilnorAlgebra_getBasis(MilnorAlgebra *algebra, uint degree);
MilnorBasisElement MilnorBasisElement_fromIndex(MilnorAlgebra *algebra, uint degree, uint index);
uint MilnorBasisElement_toIndex(MilnorAlgebra *algebra,  MilnorBasisElement b);

void MilnorAlgebra_multiply(Algebra *algebra, Vector *result, uint coeff, uint r_degree, uint r_index, uint s_degree, uint s_index);

#endif //CSTEENROD_MILNOR_H
