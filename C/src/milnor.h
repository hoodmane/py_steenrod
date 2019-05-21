//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H
#include <stdbool.h>

#include "combinatorics.h"
#include "algebra.h"

typedef char* string;

typedef struct {
    bool restricted;
    bool truncated;
    unsigned long q_part;
    unsigned long p_part_length;
    unsigned long* p_part;
} Profile;

typedef struct {
    unsigned long q_degree;
    unsigned long q_part;
    unsigned long p_degree;
    unsigned long p_length;
    unsigned long *p_part;
} MilnorBasisElement;

typedef struct {
    unsigned long length;
    MilnorBasisElement *list;
} MilnorBasisElement_list;

typedef struct {
    Algebra algebra;
    unsigned long p;
    bool generic;
    Profile profile;
    string name;
    unsigned long max_degree;
} MilnorAlgebra;


// Implemented in milnor_datatypes.c
int array_to_string(string buffer, unsigned long* A, unsigned long length);
int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b);
MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string);

// Implemented in milnor_datatypes.c
// These methods write a string to a buffer and return the length of the string written.
int milnor_element_to_string(string buffer, MilnorAlgebra * algebra, Vector * m);
int milnor_matrix_to_string(string buffer, unsigned long M[MAX_XI_TAU][MAX_XI_TAU], unsigned long rows, unsigned long cols);
int milnor_basis_element_to_key(string buffer, MilnorBasisElement *b);


MilnorAlgebra * constructMilnorAlgebra(unsigned long p, bool generic, Profile *profile);
void freeMilnorAlgebra(MilnorAlgebra *);

// Implemented in milnor.c
bool GenerateMilnorBasis(Algebra * algebra, unsigned long max_degree);
void freeMilnorBasis(MilnorAlgebra * algebra);

unsigned long GetMilnorAlgebraDimension(Algebra * algebra, unsigned long degree);
MilnorBasisElement_list GetMilnorAlgebraBasis(MilnorAlgebra * algebra, unsigned long degree);
MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, unsigned long degree, unsigned long index);
unsigned long GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);

void MilnorProduct(Algebra * algebra, Vector * result, unsigned long r_degree, unsigned long r_index, unsigned long s_degree, unsigned long s_index);

#endif //CSTEENROD_MILNOR_H
