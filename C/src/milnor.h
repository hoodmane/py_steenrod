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
int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b);
MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra *algebra, char *elt_string);

// Implemented in milnor_datatypes.c
// These methods write a string to a buffer and return the length of the string written.
int milnor_element_to_string(string buffer, MilnorAlgebra *algebra, uint degree, Vector *m);
int milnor_matrix_to_string(string buffer, uint M[MAX_XI_TAU][MAX_XI_TAU], uint rows, uint cols);
int milnor_basis_element_to_key(string buffer, MilnorBasisElement *b);


MilnorAlgebra *constructMilnorAlgebra(uint p, bool generic, Profile *profile);
void freeMilnorAlgebra(MilnorAlgebra *);

// Implemented in milnor.c
bool GenerateMilnorBasis(Algebra *algebra, uint max_degree);
void freeMilnorBasis(MilnorAlgebra *algebra);

uint GetMilnorAlgebraDimension(Algebra *algebra, uint degree);
MilnorBasisElement_list GetMilnorAlgebraBasis(MilnorAlgebra *algebra, uint degree);
MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, uint degree, uint index);
uint GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);

void MilnorProduct(Algebra *algebra, Vector *result, uint coeff, uint r_degree, uint r_index, uint s_degree, uint s_index);

#endif //CSTEENROD_MILNOR_H
