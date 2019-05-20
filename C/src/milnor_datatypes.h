//
// Created by Hood on 5/17/2019.
//

#ifndef CSTEENROD_MILNOR_DATATYPES_H
#define CSTEENROD_MILNOR_DATATYPES_H

#include <stdbool.h>
#include "khash.h"
#include "algebra.h"
typedef char* string;


typedef struct {
    string name;
    bool restricted;
    bool truncated;
    unsigned long q_part;
    unsigned long p_part_length;
    unsigned long* p_part;
} Profile;

void generate_profile_name(Profile P);
unsigned long getProfileIndex(Profile P, unsigned long index);
unsigned long getProfileExponent(Profile P, unsigned long p, unsigned long index);

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


KHASH_MAP_INIT_STR(monomial_index_map, unsigned long)

typedef struct {
    unsigned long p;
    bool generic;
    Profile profile;
    string name;
    P_part_list* P_table;
    P_part_list** P_table_by_P_length;
    unsigned long P_table_max_degree;
    Q_part_list* Q_table;
    unsigned long Q_table_max_tau;
    MilnorBasisElement_list* basis_table;
    unsigned long basis_max_degree;
    khash_t(monomial_index_map) ** basis_name_to_index_map;
} MilnorAlgebra;

void milnor_algebra_generate_name(MilnorAlgebra *A);

void milnor_algebra_generate_name(MilnorAlgebra *A);

void initializeMilnorAlgebraFields(MilnorAlgebra * A);

int array_to_string(string buffer, unsigned long* A, unsigned long length);
int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b);
MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string);
int milnor_element_to_string(string buffer, MilnorAlgebra * algebra, Vector * m);
int milnor_matrix_to_string(string buffer, unsigned long** M, unsigned long rows, unsigned long cols);

#endif //CSTEENROD_MILNOR_DATATYPES_H
