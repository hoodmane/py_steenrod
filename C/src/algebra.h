//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

#include "FpVector.h"


typedef struct Algebra {
    uint p;
    VectorInterface vectorInterface;
// Methods:
    bool (*compute_basis)(struct Algebra* this, uint degree);
    uint (*get_dimension)(struct Algebra* this, uint degree);
    void (*multiply_basis_elements)(struct Algebra* this, Vector *result, uint coeff, uint r_degree, uint r, uint s_degree, uint s);
} Algebra;

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_compute_basis(algebra, degree) (*(algebra)->compute_basis)(algebra, degree)
#define algebra_get_dimension(algebra, degree) (*(algebra)->get_dimension)(algebra, degree)
#define algebra_multiply_basis_elements(algebra, result, coeff, r_deg, r, s_deg, s) (*(algebra)->multiply_basis_elements)(algebra, result, coeff, r_deg, r, s_deg, s)





#endif //CSTEENROD_ALGEBRA_H
