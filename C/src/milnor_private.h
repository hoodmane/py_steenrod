//
// Created by Hood on 5/20/2019.
//

#ifndef C_MILNOR_PRIVATE_H
#define C_MILNOR_PRIVATE_H

#include "khash.h"

#include "milnor.h"

typedef struct {
    uint degree;
    uint length;
    uint *p_part;
} P_part;

typedef struct {
    uint length;
    P_part *list;
} P_part_list;

typedef struct {
    uint degree;
    uint bit_string;
} Q_part;

typedef struct {
    uint length;
    Q_part *list;
} Q_part_list;

KHASH_MAP_INIT_STR(monomial_index_map, uint)

typedef struct {
    MilnorAlgebra public_algebra;
    P_part_list* P_table;
    P_part_list** P_table_by_P_length;
    uint P_table_max_degree;
    Q_part_list* Q_table;
    uint Q_table_max_tau;
    MilnorBasisElement_list* basis_table;
    khash_t(monomial_index_map) ** basis_element_to_index_map;
} MilnorAlgebraInternal;

void milnor_algebra_generate_name(MilnorAlgebra *A);
uint Profile_getExponent(Profile P, uint p, uint index);

// Private functions
void generateMilnorBasisPpartTable(MilnorAlgebraInternal * algebra, uint n);
void freeMilnorBasisPpartTable(MilnorAlgebraInternal * algebra);
void generateMilnorBasisQpartTable(MilnorAlgebraInternal * algebra, uint n );
void freeMilnorBasisQPartTable(MilnorAlgebraInternal * algebra);

void GenerateMilnorBasis2(MilnorAlgebraInternal * algebra, uint old_n, uint n);
void GenerateMilnorBasisGeneric(MilnorAlgebraInternal * algebra, uint old_n, uint n);

void MilnorProductEven(MilnorAlgebraInternal * algebra, uint * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
void MilnorProductFullQpart(MilnorAlgebraInternal * algebra, uint * result, MilnorBasisElement m1, uint f);

void initialize_milnor_matrix(uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);
bool step_milnor_matrix(uint  p, uint M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);

#endif //C_MILNOR_PRIVATE_H
