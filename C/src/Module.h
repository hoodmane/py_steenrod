#ifndef CSTEENROD_MODULE_H
#define CSTEENROD_MODULE_H

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "Algebra.h"

#define MODULE_TYPE_FREE 0
#define MODULE_TYPE_FINITE_DIMENSIONAL 1
#define MODULE_TYPE_FINITELY_PRESENTED 2

typedef struct Module {
    uint p;
    Algebra *algebra;    
    uint type;
    int min_degree;
    int max_degree;        
// Methods:
    bool (*computeBasis)(struct Module *this, int degree);
    uint (*getDimension)(struct Module *this, int degree);
    void (*actOnBasis)(struct Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index);
} Module;

#define Module_computeBasis(module, degree) ((module)->computeBasis)(module, degree)
#define Module_getDimension(module, degree) ((module)->getDimension)(module, degree)
#define Module_actOnBasis(module, result, coeff, op_deg, op, r_deg, r) ((module)->actOnBasis)(module, result, coeff, op_deg, op, r_deg, r)

// For javascript
uint Module_getDimension_function(Module *module, int degree);




#endif //CSTEENROD_MODULE_H