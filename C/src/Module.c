//
// Created by Hood on 5/20/2019.
//

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "FpVector.h"
#include "Algebra.h"
#include "Module.h"

// For javascript
bool Module_computeBasis_function(Module *this, int degree){
    return Module_computeBasis(this, degree);
}

uint Module_getDimension_function(Module *module, int degree){
    uint result = module->getDimension(module, degree);
    return result;
}

void Module_actOnBasis_function(Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index){
    Module_actOnBasis(this, result, coeff, op_degree, op_index, mod_degree, mod_index);
}



/*
#include "milnor.h"
int main(){
    initializePrime(2);
    Algebra *algebra = (Algebra*)constructMilnorAlgebra(2, false, NULL);
    GenerateMilnorBasis(algebra, 50);
    // FreeModule *F = FreeModule_construct(algebra, 2, 20);
    // F->number_of_generators = 2;
    // F->number_of_generators_in_degree[0] = 1;
    // F->number_of_generators_in_degree[1] = 1;
    // Module *Fm = (Module*) F;
    // FreeModule_ConstructBlockOffsetTable(F, 0);
    // FreeModule_ConstructBlockOffsetTable(F, 1);
    // FreeModule_ConstructBlockOffsetTable(F, 2);
    // FreeModule_ConstructBlockOffsetTable(F, 3);
    // Vector *v = constructVector2(2, 3, 0);
    // FreeModuleInternal *Fi = (FreeModuleInternal*)F;
    // printf("gen_indices: %d, ", Fi->basis_element_to_opgen_table[1][0].generator_degree);
    // printf("%d, \n", Fi->basis_element_to_opgen_table[1][1].generator_degree);
    // FreeModule_actOnBasis(Fm, v, 1, 3, 0, 0, 0);
    // FreeModule_actOnBasis(Fm, v, 1, 2, 0, 1, 0);
    // FreeModule_actOnBasis(Fm, v, 1, 2, 0, 1, 1);
    // printVector(v);

    FreeModule *F1 = FreeModule_construct(algebra, 4, 20);
    FreeModule *F0 = FreeModule_construct(algebra, 0, 20);
    F0->number_of_generators = 1;
    F0->number_of_generators_in_degree[0] = 1;
    for(uint i=0; i<10; i++){
        FreeModule_ConstructBlockOffsetTable(F0, i);
    }
    F1->number_of_generators_in_degree[1] = 1;
    F1->number_of_generators_in_degree[2] = 1;
    F1->number_of_generators_in_degree[4] = 1;
    FreeModuleHomomorphism *f = FreeModuleHomomorphism_construct(F1, (Module*)F0, 20);
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f, 1, 1);
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f, 2, 1);
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f, 4, 1);
    Vector *output1 = constructVector2(2, 1, 0);
    uint array[1];
    array[0] = 1;
    packVector(output1, array);
    Vector *output2 = constructVector2(2, 1, 0);
    packVector(output2, array);
    uint array4[2];
    array4[0] = 1;
    array4[1] = 0;
    Vector *output4 = constructVector2(2, 2, 0);
    packVector(output4, array4);
    FreeModuleHomomorphism_setOutput(f, 1, 0, output1);
    FreeModuleHomomorphism_setOutput(f, 2, 0, output2);
    FreeModuleHomomorphism_setOutput(f, 4, 0, output4);
    uint degree = 5;
    Matrix *result = constructMatrix2(2, 
        module_getDimension((Module*)F1, degree), 
        module_getDimension((Module*)F0, degree)
    );
    FreeModuleHomomorphism_getMatrix(f, result, degree);    
    printf("result(%d): \n", degree); printMatrix(result);
}
*/