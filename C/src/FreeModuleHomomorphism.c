#include <assert.h>

#include "FreeModule.h"
#include "FreeModuleHomomorphism.h"

// max_degree is one larger than the highest degree in which you are allowed to access this module.
FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule *source, Module *target, int degree_shift, int max_degree){
    uint num_degrees = max_degree - source->module.min_degree;
    FreeModuleHomomorphism *f = malloc(
        sizeof(FreeModuleHomomorphism) 
        + num_degrees * sizeof(Vector**) // outputs
        + num_degrees * sizeof(Matrix*)  // coimage_to_image_iso
        + num_degrees * sizeof(Subspace*)  // kernel
    );
    f->source = source;
    f->target = target;
    f->max_degree = max_degree;
    f->max_computed_degree = source->module.min_degree;
    f->degree_shift = degree_shift;
    f->outputs = (Vector***)(f + 1);
    f->coimage_to_image_isomorphism = (Matrix**)(f->outputs + num_degrees);
    f->kernel = (Subspace**)(f->coimage_to_image_isomorphism + num_degrees);
    memset(f+1, 0, num_degrees * (sizeof(Vector**) + sizeof(Matrix*) + sizeof(Subspace*)));
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
        Subspace_free(f->kernel[i]);
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
    uint dimension = Module_getDimension(f->target, degree + f->degree_shift);
    uint vector_size = Vector_getSize(p, dimension, 0);
    f->outputs[shifted_degree] = (Vector**)malloc(
        num_gens * sizeof(Vector*) 
        + num_gens * vector_size
    );
    Vector **vector_ptr_ptr = f->outputs[shifted_degree];
    char *vector_memory_ptr = (char*)(vector_ptr_ptr + num_gens);
    for(uint i = 0; i < num_gens; i++){
        f->outputs[shifted_degree][i] = Vector_initialize(p, &vector_memory_ptr, dimension, 0);
        vector_ptr_ptr ++;
    }
    assert(vector_ptr_ptr == f->outputs[shifted_degree] + num_gens);
    assert(vector_memory_ptr == (char*)(vector_ptr_ptr) + num_gens * vector_size);
}

Vector *FreeModuleHomomorphism_getOutput(FreeModuleHomomorphism *f, int generator_degree, uint generator_index){
    assert(generator_degree - f->source->module.min_degree >= 0);
    assert(generator_degree < f->max_computed_degree);
    assert(generator_index < f->source->number_of_generators_in_degree[generator_degree]);
    return f->outputs[generator_degree - f->source->module.min_degree][generator_index];
}

void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, int generator_degree, uint generator_index, Vector *output){
    assert(output->dimension == Module_getDimension(f->target, generator_degree + f->degree_shift));
    assert(output->offset == 0);
    assert(generator_index < FreeModule_getNumberOfGensInDegree(f->source, generator_degree));
    Vector_assign(f->outputs[generator_degree - f->source->module.min_degree][generator_index], output);
}

void FreeModuleHomomorphism_addGeneratorsFromMatrixRows(FreeModuleHomomorphism *f, uint degree, Matrix* matrix, uint first_new_row, uint new_generators){
    uint p = matrix->p;
    uint dimension = Module_getDimension(f->target, degree);
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(f, degree, new_generators);
    char slice_memory[Vector_getSize(p, 0, 0)]; 
    char *slice_ptr = slice_memory;
    Vector *slice = Vector_initialize(p, &slice_ptr, 0, 0);
    for(uint i = 0; i < new_generators; i++){
        Vector_slice(slice, matrix->vectors[first_new_row + i], 0, dimension);
        FreeModuleHomomorphism_setOutput(f, degree, i, slice);
    }
}

// Primarily for Javascript (so we can avoid indexing struct fields).
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
    assert(input_index < Module_getDimension((Module*)f->source, input_degree));
    assert(Module_getDimension(f->target, input_degree + f->degree_shift) == result->dimension);
    FreeModuleOperationGeneratorPair operation_generator = 
        FreeModule_indexToOpGen(f->source, input_degree, input_index);
    int operation_degree = operation_generator.operation_degree;
    uint operation_index = operation_generator.operation_index;
    int generator_degree = operation_generator.generator_degree;
    uint generator_index = operation_generator.generator_index;
    Vector *output_on_generator = FreeModuleHomomorphism_getOutput(f, generator_degree, generator_index);
    Module_actOnElement(f->target, result, coeff, operation_degree, operation_index, generator_degree + f->degree_shift, output_on_generator);
}

void FreeModuleHomomorphism_apply(FreeModuleHomomorphism *f, Vector *result, uint coeff, int input_degree, Vector *input){
    assert(input->dimension == Module_getDimension((Module*)f->source, input_degree));
    for(        
        VectorIterator it = Vector_getIterator(input);
        it.has_more; 
        it = Vector_stepIterator(it)
    ){
        if(it.value!=0){
            FreeModuleHomomorphism_applyToBasisElement(f, result, coeff*it.value, input_degree, it.index);
        }
    }
}

// result should be big enough to hold output (how big is that?)
// Well it should have dim(target) columns and dim(source) rows. I guess this is 
// the transpose of the usual convention.
void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, int degree){
    assert(degree < f->max_degree);
    assert(Module_getDimension(&f->source->module, degree) <= result->rows);
    assert(Module_getDimension(f->target, degree + f->degree_shift) <= result->columns);
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
                Module_actOnElement(f->target, result->vectors[i], 1, op_deg, op_idx, gen_deg + f->degree_shift, output_on_generator);
                i++;
            }
        }
    }
}

void FreeModuleHomomorphism_computeKernel(FreeModuleHomomorphism *f, Matrix *full_matrix, int *pivots, uint degree){
    uint p = full_matrix->p;
    uint source_dimension = Module_getDimension((Module*)f->source, degree);
    uint target_dimension = Module_getDimension(f->target, degree);
    uint first_source_index = Vector_getPaddedDimension(p, target_dimension, 0);        
    // For the first stage we just want the part of size source_dimension x (padded_target_dimension + source_dimension).
    // Slice matrix out of full_matrix.
    char slice_matrix_memory[Matrix_getSliceSize(full_matrix->p, source_dimension)];
    Matrix *matrix = Matrix_slice(full_matrix, slice_matrix_memory, 0, source_dimension, 0, first_source_index + source_dimension);
    FreeModuleHomomorphism_getMatrix(f, matrix, degree);

    // Write the identity matrix into the right block
    for(uint i = 0; i < source_dimension; i++){
        Vector_setEntry(matrix->vectors[i], first_source_index + i, 1);
    }

    // Row reduce
    rowReduce(matrix, pivots, 0, 0);//target_dimension, padded_target_dimension);
    // Make sure to permute the rows of the full_matrix so they are consistent with the slice.
    uint permutation[matrix->rows];
    Matrix_getRowPermutation(matrix, permutation);
    Matrix_applyRowPermutation(full_matrix, permutation, matrix->rows); 

    Subspace *kernel = Matrix_computeKernel(matrix, pivots, first_source_index);
    f->kernel[degree - f->target->min_degree] = kernel;
}


FreeModule *FreeModuleHomomorphism_getSource(FreeModuleHomomorphism *f){
    return f->source;
}

Module *FreeModuleHomomorphism_getTarget(FreeModuleHomomorphism *f){
    return f->target;
}