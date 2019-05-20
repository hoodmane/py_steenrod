//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

typedef struct Algebra {
    unsigned long p;
// Methods:
    bool (*compute_basis)(struct Algebra* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Algebra* this, unsigned long degree);
    int (*multiply_basis_elements)(struct Algebra* this, unsigned long degree, unsigned long index1, unsigned long index2);
} Algebra;


typedef struct Module {
    unsigned long p;
    Algebra * algebra;
// Methods:
    bool (*compute_basis)(struct Module* this, unsigned long degree);
    unsigned long (*get_basis_dimension)(struct Module* this, unsigned long degree);
    int (*act_on_basis)(struct Module* this, unsigned long degree, unsigned long alg_index, unsigned long mod_index);
} Module;

#endif //CSTEENROD_ALGEBRA_H
