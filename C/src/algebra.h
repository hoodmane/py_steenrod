//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

#include "FpVector.h"


typedef struct Algebra {
    uint p;
// Methods:
    bool (*computeBasis)(struct Algebra* this, uint degree);
    uint (*getDimension)(struct Algebra* this, uint degree);
    void (*multiplyBasisElements)(struct Algebra* this, Vector *result, uint coeff, uint r_degree, uint r, uint s_degree, uint s);
} Algebra;

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_computeBasis(algebra, degree) (*(algebra)->computeBasis)(algebra, degree)
#define algebra_getDimension(algebra, degree) (*(algebra)->getDimension)(algebra, degree)
#define algebra_multiplyBasisElements(algebra, result, coeff, r_deg, r, s_deg, s) (*(algebra)->multiplyBasisElements)(algebra, result, coeff, r_deg, r, s_deg, s)





#endif //CSTEENROD_ALGEBRA_H
