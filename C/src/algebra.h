//
// Created by Hood on 5/20/2019.
//

#ifndef CSTEENROD_ALGEBRA_H
#define CSTEENROD_ALGEBRA_H
#include <stdbool.h>

#include "FpVector.h"

// List of structlines to indicate.
// TODO: generalize to higher filtration elements?
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
#define Algebra_computeBasis(algebra, degree) (*(algebra)->computeBasis)(algebra, degree)
#define Algebra_getDimension(algebra, degree, excess) (*(algebra)->getDimension)(algebra, degree, excess)
#define Algebra_multiplyBasisElements(algebra, result, coeff, r_deg, r, s_deg, s, excess) (*(algebra)->multiplyBasisElements)(algebra, result, coeff, r_deg, r, s_deg, s, excess)
#define Algebra_basisElementToString(algebra, result, degree, idx) (*(algebra)->basisElementToString)(algebra, result, degree, idx)

// These are for calling from python / javascript
void Algebra_computeBasis_function(Algebra *algebra, int degree);
uint Algebra_getDimension_function(Algebra *algebra, int degree, int excess);
void Algebra_multiplyBasisElements_function(Algebra *algebra, Vector *result, uint coeff, int r_deg, uint r_idx, int s_deg, uint s_idx, int excess);
uint Algebra_basisElementToString_function(Algebra *algebra, char *result, int degree, uint idx);



#endif //CSTEENROD_ALGEBRA_H
