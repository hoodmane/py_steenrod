//
// Created by Hood on 5/20/2019.
//

#ifndef C_MILNOR_PRIVATE_H
#define C_MILNOR_PRIVATE_H

#include "milnor.h"
#include "khash.h"

typedef struct {
    unsigned long degree;
    unsigned long length;
    unsigned long *p_part;
} P_part;

typedef struct {
    unsigned long length;
    P_part *list;
} P_part_list;

typedef struct {
    unsigned long degree;
    unsigned long bit_string;
} Q_part;

typedef struct {
    unsigned long length;
    Q_part *list;
} Q_part_list;

KHASH_MAP_INIT_STR(monomial_index_map, unsigned long)

typedef struct {
    MilnorAlgebra public_algebra;
    P_part_list* P_table;
    P_part_list** P_table_by_P_length;
    unsigned long P_table_max_degree;
    Q_part_list* Q_table;
    unsigned long Q_table_max_tau;
    MilnorBasisElement_list* basis_table;
    khash_t(monomial_index_map) ** basis_element_to_index_map;
} MilnorAlgebraInternal;

void milnor_algebra_generate_name(MilnorAlgebra *A);
unsigned long getProfileExponent(Profile P, unsigned long p, unsigned long index);

// Private functions
void generateMilnorBasisPpartTable(MilnorAlgebraInternal * algebra, unsigned long n);
void freeMilnorBasisPpartTable(MilnorAlgebraInternal * algebra);
void generateMilnorBasisQpartTable(MilnorAlgebraInternal * algebra, unsigned long n );
void freeMilnorBasisQPartTable(MilnorAlgebraInternal * algebra);

void GenerateMilnorBasis2(MilnorAlgebraInternal * algebra, unsigned long old_n, unsigned long n);
void GenerateMilnorBasisGeneric(MilnorAlgebraInternal * algebra, unsigned long old_n, unsigned long n);

void MilnorProductEven(MilnorAlgebraInternal * algebra, Vector * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
void MilnorProductFullQpart(MilnorAlgebraInternal * algebra, Vector * result, MilnorBasisElement m1, unsigned long f);

void initialize_milnor_matrix(unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);
bool step_milnor_matrix(unsigned long  p, unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);

#endif //C_MILNOR_PRIVATE_H
