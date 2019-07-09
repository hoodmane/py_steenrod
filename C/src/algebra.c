#include "algebra.h"

void Algebra_computeBasis_function(Algebra *algebra, int degree){
    Algebra_computeBasis(algebra, degree);
}

uint Algebra_getDimension_function(Algebra *algebra, int degree, int excess){
    return Algebra_getDimension(algebra, degree, excess);
}

void Algebra_multiplyBasisElements_function(Algebra *algebra, Vector *result, uint coeff, int r_deg, uint r_idx, int s_deg, uint s_idx, int excess){
    Algebra_multiplyBasisElements(algebra, result, coeff, r_deg, r_idx, s_deg, s_idx, excess);
}

uint Algebra_basisElementToString_function(Algebra *algebra, char *result, int degree, uint idx){
    return Algebra_basisElementToString(algebra, result, degree, idx);
}