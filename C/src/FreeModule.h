#ifndef CSTEENROD_FREE_MODULE_H
#define CSTEENROD_FREE_MODULE_H
#include "Module.h"

typedef struct {
    Module module;
    uint *number_of_generators_in_degree;
} FreeModule;

typedef struct {
    int operation_degree;
    uint operation_index;
    int generator_degree;
    uint generator_index;
} FreeModuleOperationGeneratorPair;

FreeModule *FreeModule_construct(Algebra *algebra, int min_degree, int max_degree);
void FreeModule_free(FreeModule *module);

bool FreeModule_computeBasis(Module *this, int degree);
uint FreeModule_getNumberOfGensInDegree(FreeModule *this, int degree);
uint FreeModule_getDimension(Module *this, int degree);
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, 
    int op_degree, uint op_index, int mod_degree, uint mod_idx);

void FreeModule_ConstructBlockOffsetTable(FreeModule *M, int degree);
void FreeModule_addGenerators(FreeModule *this, int degree, uint new_generators);
uint FreeModule_operationGeneratorToIndex(FreeModule *this, int op_deg, uint op_idx, int gen_deg, uint gen_idx);
FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, int degree, uint index);

uint FreeModule_element_toJSONString(char *result, FreeModule *this, int degree, Vector *element);


#endif //CSTEENROD_FREE_MODULE_H