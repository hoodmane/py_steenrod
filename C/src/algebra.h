//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

#include "FpVector.h"

typedef struct {
    uint length;
    int *degrees;
    uint *indices;
} FiltrationOneProductList;

typedef struct Algebra {
    uint p;
    uint max_degree; 
    char *name;
    FiltrationOneProductList *product_list; // This determines which indecomposibles have lines drawn for them.
// Methods:
    void (*computeBasis)(struct Algebra *this, int degree);
    uint (*getDimension)(struct Algebra *this, int degree, int excess);
    void (*multiplyBasisElements)(struct Algebra *this, Vector *result, uint coeff, int r_degree, uint r_idx, int s_degree, uint s_idx, int excess);
    uint (*basisElementToString)(struct Algebra *this, char *result, int degree, uint idx);
} Algebra;

// Careful with these macros: could cause multiple evaluation of algebra / module.
#define algebra_computeBasis(algebra, degree) (*(algebra)->computeBasis)(algebra, degree)
#define algebra_getDimension(algebra, degree, excess) (*(algebra)->getDimension)(algebra, degree, excess)
#define algebra_multiplyBasisElements(algebra, result, coeff, r_deg, r, s_deg, s, excess) (*(algebra)->multiplyBasisElements)(algebra, result, coeff, r_deg, r, s_deg, s, excess)
#define algebra_basisElementToString(algebra, result, degree, idx) (*(algebra)->basisElementToString)(algebra, result, degree, idx)

// These are for calling from python / javascript
void algebra_computeBasis_function(Algebra *algebra, int degree);
uint algebra_getDimension_function(Algebra *algebra, int degree, int excess);
void algebra_multiplyBasisElements_function(Algebra *algebra, Vector *result, uint coeff, int r_deg, uint r_idx, int s_deg, uint s_idx, int excess);
uint algebra_basisElementToString_function(Algebra *algebra, char *result, int degree, uint idx);



#endif //CSTEENROD_ALGEBRA_H
