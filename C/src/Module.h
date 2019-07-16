#ifndef CSTEENROD_MODULE_H
#define CSTEENROD_MODULE_H

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "Algebra.h"

// TODO: Maybe use these?
#define MODULE_TYPE_FREE 0
#define MODULE_TYPE_FINITE_DIMENSIONAL 1
#define MODULE_TYPE_FINITELY_PRESENTED 2

typedef struct Module {
    uint p;
    Algebra *algebra;
    char *name;
    uint type;
    int min_degree;
    int max_degree; // Rename to max_allocated_degree?
    int max_computed_degree;
// Methods:
    bool (*computeBasis)(struct Module *this, int degree);
    uint (*getDimension)(struct Module *this, int degree);
    void (*actOnBasis)(struct Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index);
} Module;

#define Module_computeBasis(module, degree) ((module)->computeBasis)(module, degree)
#define Module_getDimension(module, degree) ((module)->getDimension)(module, degree)
#define Module_actOnBasis(module, result, coeff, op_deg, op, r_deg, r) ((module)->actOnBasis)(module, result, coeff, op_deg, op, r_deg, r)

void Module_actOnElement(Module *this, Vector *result, uint coeff, int op_deg, uint op_idx, int module_degree, Vector *module_element);

// For javascript
bool Module_computeBasis_function(Module *this, int degree);
uint Module_getDimension_function(Module *this, int degree);
void Module_actOnBasis_function(Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index);




#endif //CSTEENROD_MODULE_H