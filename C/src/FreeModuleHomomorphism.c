#include "FreeModule.h"
#include "FreeModuleHomomorphism.h"
#include <assert.h>

// max_degree is one larger than the highest degree in which you are allowed to access this module.
FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule *source, Module *target, int max_degree){
    uint num_degrees = max_degree - source->module.min_degree;
    FreeModuleHomomorphism *f = malloc(
        sizeof(FreeModuleHomomorphism) 
        + num_degrees * sizeof(Vector**) // outputs
        + num_degrees * sizeof(Matrix*)  // coimage_to_image_iso
        + num_degrees * sizeof(Kernel*)  // kernel
    );
    f->source = source;
    f->target = target;
    f->max_degree = max_degree;
    f->max_computed_degree = source->module.min_degree;
    f->outputs = (Vector***)(f + 1);
    f->coimage_to_image_isomorphism = (Matrix**)(f->outputs + num_degrees);
    f->kernel = (Kernel**)(f->coimage_to_image_isomorphism + num_degrees);
    memset(f+1, 0, num_degrees * (sizeof(Vector**) + sizeof(Matrix*) + sizeof(Kernel*)));
    assert(f->coimage_to_image_isomorphism[0] == NULL);
    assert(f->kernel[0] == NULL);
    return f;
}

void FreeModuleHomomorphism_free(FreeModuleHomomorphism *f){
    if(f == NULL){
        return;
    }
    for(int i = 0; i < f->max_degree - f->source->module.min_degree; i++){
        Matrix_free(f->coimage_to_image_isomorphism[i]);
        Kernel_free(f->kernel[i]);
        free(f->outputs[i]);
    }
    free(f);
}

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, int degree, uint num_gens){
    assert(degree < f->max_degree);
    assert(f->max_computed_degree <= degree);
    int shifted_degree = degree - f->source->module.min_degree;
    assert(shifted_degree >= 0);
    f->max_computed_degree = degree + 1;
    uint p = f->source->module.p;
    uint dimension = Module_getDimension(f->target, degree);
    uint vector_size = Vector_getSize(p, dimension, 0);
    f->outputs[shifted_degree] = (Vector**)malloc(
        num_gens * sizeof(Vector*) 
        + num_gens * Vector_getContainerSize(p)
        + num_gens * vector_size
    );
    Vector **vector_ptr_ptr = f->outputs[shifted_degree];
    char *vector_container_ptr = (char*)(vector_ptr_ptr + num_gens);
    char *vector_memory_ptr = vector_container_ptr + num_gens * Vector_getContainerSize(p);
    for(uint i = 0; i < num_gens; i++){
        f->outputs[shifted_degree][i] = Vector_initialize(p, vector_container_ptr, vector_memory_ptr, dimension, 0);
        vector_ptr_ptr ++;
        vector_container_ptr += Vector_getContainerSize(p);
        vector_memory_ptr += vector_size;
    }
    assert(vector_ptr_ptr == f->outputs[shifted_degree] + num_gens);
    assert(vector_container_ptr == (char*)(vector_ptr_ptr) + num_gens * Vector_getContainerSize(p));
    assert(vector_memory_ptr == vector_container_ptr + num_gens * vector_size);

}

Vector *FreeModuleHomomorphism_getOutput(FreeModuleHomomorphism *f, int generator_degree, uint generator_index){
    assert(generator_degree - f->source->module.min_degree >= 0);
    assert(generator_degree < f->max_computed_degree);
    assert(generator_index < Module_getDimension((Module*)f->source, generator_degree));
    return f->outputs[generator_degree - f->source->module.min_degree][generator_index];
}

void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, int generator_degree, uint generator_index, Vector *output){
    assert(output->dimension == Module_getDimension(f->target, generator_degree));
    assert(output->offset == 0);
    assert(generator_index < FreeModule_getNumberOfGensInDegree(f->source, generator_degree));
    Vector_assign(f->outputs[generator_degree - f->source->module.min_degree][generator_index], output);
}

// Primarily for Javascript (so we can avoid actually indexing struct fields).
void FreeModuleHomomorphism_applyToGenerator(FreeModuleHomomorphism *f, Vector *result, uint coeff, int generator_degree, uint generator_index){
    Vector *output_on_generator = FreeModuleHomomorphism_getOutput(f, generator_degree, generator_index);
    Vector_add(result, output_on_generator, coeff);
}

// Run FreeModule_ConstructBlockOffsetTable(source, degree) before using this on an input in that degree
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, int input_degree, uint input_index){
    int shifted_input_degree = input_degree - f->source->module.min_degree;
    assert(shifted_input_degree > 0);
    assert(input_degree < f->max_degree);
    // assert(((FreeModuleInternal*)f->source)->basis_element_to_opgen_table[shifted_input_degree] != NULL);
    assert(input_index < FreeModule_getDimension((Module*)f->source, input_degree));
    FreeModuleOperationGeneratorPair operation_generator = 
        FreeModule_indexToOpGen(f->source, input_degree, input_index);
    int operation_degree = operation_generator.operation_degree;
    uint operation_index = operation_generator.operation_index;
    int generator_degree = operation_generator.generator_degree;
    uint generator_index = operation_generator.generator_index;
    Vector *output_on_generator = FreeModuleHomomorphism_getOutput(f, generator_degree, generator_index);
    for(
        VectorIterator it = Vector_getIterator(output_on_generator);
        it.has_more; 
        it = Vector_stepIterator(it)
    ){
        if(it.value != 0){
            uint c = modPLookup( f->source->module.p, it.value*coeff);
            Module_actOnBasis(f->target, result, c, operation_degree, operation_index, generator_degree, it.index);
        }
    }
}

// result should be big enough to hold output (how big is that?)
// Well it should have dim(target) columns and dim(source) rows. I guess this is 
// the transpose of the usual convention.
void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, int degree){
    assert(degree < f->max_degree);
    assert(Module_getDimension(&f->source->module, degree) <= result->rows);
    assert(Module_getDimension(f->target, degree) <= result->columns);
    // The shorter implementation if we do FreeModuleConstructBlockOffsetTable first.
    // Maybe we ought to do that...
    // for(int i = 0; i < Module_getDimension(&f->source->module, degree); i++){
    //     FreeModuleHomomorphism_applyToBasisElement(f, result->matrix[i], 1, degree, i);
    // }    
    FreeModule *source = f->source;
    Algebra *algebra = source->module.algebra;
    // 
    uint i = 0;
    for(int gen_deg = f->source->module.min_degree; gen_deg <= degree; gen_deg++){
        int op_deg = degree - gen_deg;
        uint num_ops = Algebra_getDimension(algebra, op_deg, gen_deg);
        for(uint gen_idx = 0; gen_idx < FreeModule_getNumberOfGensInDegree(f->source, gen_deg); gen_idx++){
            for(uint op_idx = 0; op_idx < num_ops; op_idx++){
                Vector *output_on_generator = FreeModuleHomomorphism_getOutput(f, gen_deg, gen_idx);
                for(
                    VectorIterator it = Vector_getIterator(output_on_generator);
                    it.has_more; 
                    it = Vector_stepIterator(it)
                ){
                    if(it.value != 0){
                        // our element of our source is op * gen. It maps to op * (f(gen)).
                        Module_actOnBasis(f->target, result->matrix[i], it.value, op_deg, op_idx, gen_deg, it.index);
                    }
                }
                i++;                
            }
        }
    }
}

FreeModule *FreeModuleHomomorphism_getSource(FreeModuleHomomorphism *f){
    return f->source;
}

Module *FreeModuleHomomorphism_getTarget(FreeModuleHomomorphism *f){
    return f->target;
}