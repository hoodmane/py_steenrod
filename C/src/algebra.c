#include "algebra.h"

void algebra_computeBasis_function(Algebra *algebra, int degree){
    algebra_computeBasis(algebra, degree);
}

uint algebra_getDimension_function(Algebra *algebra, int degree, int excess){
    return algebra_getDimension(algebra, degree, excess);
}

void algebra_multiplyBasisElements_function(Algebra *algebra, Vector *result, uint coeff, int r_deg, uint r_idx, int s_deg, uint s_idx, int excess){
    algebra_multiplyBasisElements(algebra, result, coeff, r_deg, r_idx, s_deg, s_idx, excess);
}

uint algebra_basisElementToString_function(Algebra *algebra, char *result, int degree, uint idx){
    return algebra_basisElementToString(algebra, result, degree, idx);
}