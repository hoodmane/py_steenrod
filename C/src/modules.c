//
// Created by Hood on 5/20/2019.
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "modules.h"
#include "FpVector.h"
#include "algebra.h"


// The allocator is horrendous so we're going to separate it out.
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint max_basis_degree, uint *graded_dimension);

FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, int min_degree, int max_basis_degree, uint *graded_dimension){
    FiniteDimensionalModule *result = FiniteDimensionalModule_allocate(algebra, max_basis_degree - min_degree, graded_dimension);
    result->module.p = algebra->p;
    result->module.algebra = algebra;
    result->module.computeBasis = FiniteDimensionalModule_computeBasis;
    result->module.getDimension = FiniteDimensionalModule_getDimension;
    result->module.actOnBasis = FiniteDimensionalModule_actOnBasis;
    result->module.min_degree = min_degree;
    result->module.max_degree = -1; // There is no "max degree" for a finite module -- we've computed it through an infinite range.
    result->max_basis_degree = max_basis_degree;
    memcpy(result->graded_dimension, graded_dimension, (max_basis_degree - min_degree) * sizeof(uint));
    return result;
}

void FiniteDimensionalModule_free(FiniteDimensionalModule *module){
    free(module);
}

// This is the grossest allocator.
// The most important part of the FD module is the 5d ragged array used to fetch actions.
// It takes a fair bit of effort to lay it out in memory...
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint max_basis_degree, uint *graded_dimension){
    uint p = algebra->p;
    uint vector_container_size = Vector_getContainerSize(p);
    // Count number of triples (x, y, op) with |x| + |op| = |y|.
    // The amount of memory we need to allocate is:
    // # of input_degrees  * sizeof(***Vector)
    // + # of nonempty input degrees * # of output degrees * sizeof(**Vector)
    // + Sum over (nonempty in_deg < nonempty out_deg) of (
    //              # of operations in (out_deg - in_deg) * sizeof(*Vector)
    //              # of operations in (out_deg - in_deg) * # of gens in degree in_degree * sizeof(Vector)
    //              # of operations in (out_deg - in_deg) * # of gens in degree in_degree * # of gens in degree out_degree * sizeof(uint)
    // )
    size_t action_matrix_size_1 = 0, action_matrix_size_2 = 0,
            action_matrix_size_3 = 0, action_matrix_size_4 = 0,
            action_matrix_size_5 = 0;
    uint number_of_nonempty_degrees = 0;
    for(int degree = 0; degree < max_basis_degree; degree++) {
        if(graded_dimension[degree] != 0){
            number_of_nonempty_degrees ++;
        }
    }
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> (out_index) -> value
    //  ****    -> ***       -> **Vector   -> *Vector    -> Vector -> uint
    action_matrix_size_1 += max_basis_degree * sizeof(Vector***);
    action_matrix_size_2 += number_of_nonempty_degrees * max_basis_degree * sizeof(Vector**);
    for(uint input_degree = 0; input_degree < max_basis_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            continue;
        }
        for(uint output_degree = input_degree + 1; output_degree < max_basis_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                continue;
            }
            uint number_of_operations = algebra_getDimension(algebra, output_degree - input_degree, input_degree);
            action_matrix_size_3 += sizeof(Vector**) * number_of_operations;
            action_matrix_size_4 += sizeof(Vector*) * number_of_operations * graded_dimension[input_degree];
            uint vector_size = Vector_getSize(p, graded_dimension[output_degree], 0);
            action_matrix_size_5 += (vector_container_size + vector_size) * number_of_operations * graded_dimension[input_degree];
        }
    }

    action_matrix_size_2 += action_matrix_size_1;
    action_matrix_size_3 += action_matrix_size_2;
    action_matrix_size_4 += action_matrix_size_3;
    action_matrix_size_5 += action_matrix_size_4;

    FiniteDimensionalModule *result = malloc(
            sizeof(FiniteDimensionalModule)
            + max_basis_degree * sizeof(uint) // graded_dimension
            + action_matrix_size_5
    );
    result->graded_dimension = (uint *) (result + 1);
    char *top_of_action_table = (char *) (result + 1) + max_basis_degree *sizeof(uint);
    Vector *****current_ptr_1 = (Vector *****) top_of_action_table;
    Vector ****current_ptr_2 = (Vector ****) (top_of_action_table + action_matrix_size_1);
    Vector ***current_ptr_3 = (Vector ***) (top_of_action_table + action_matrix_size_2);
    Vector **current_ptr_4 = (Vector **) (top_of_action_table + action_matrix_size_3);
    char *current_ptr_5 = (char *) (top_of_action_table + action_matrix_size_4);
    result->actions = current_ptr_1; 
    memset(top_of_action_table, 0, action_matrix_size_4);
    for(uint input_degree = 0; input_degree < max_basis_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            current_ptr_1 ++;
            continue;
        }
        *current_ptr_1 = current_ptr_2;
        current_ptr_2 += input_degree + 1;
        for(uint output_degree = input_degree + 1; output_degree < max_basis_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                current_ptr_2 ++;
                continue;
            }
            *current_ptr_2 = current_ptr_3;
            uint vector_size = Vector_getSize(p, graded_dimension[output_degree], 0);
            uint vector_total_size = vector_size + vector_container_size;
            uint number_of_operations = algebra_getDimension(algebra, output_degree - input_degree, input_degree);
            for(uint operation_idx = 0; operation_idx < number_of_operations; operation_idx ++){
                *current_ptr_3 = current_ptr_4;
                for(uint input_idx = 0; input_idx < graded_dimension[input_degree]; input_idx ++ ){
                    *current_ptr_4 = Vector_initialize(p, current_ptr_5, current_ptr_5 + vector_container_size, graded_dimension[output_degree], 0);
                    current_ptr_4 ++;
                    current_ptr_5 += vector_total_size;
                }
                current_ptr_3 ++;
            }
            current_ptr_2 ++;
        }
        current_ptr_1 ++;
    }
    assert((uint64)current_ptr_1 == (uint64)(top_of_action_table + action_matrix_size_1));
    assert((uint64)current_ptr_2 == (uint64)(top_of_action_table + action_matrix_size_2));
    assert((uint64)current_ptr_3 == (uint64)(top_of_action_table + action_matrix_size_3));
    assert((uint64)current_ptr_4 == (uint64)(top_of_action_table + action_matrix_size_4));
    assert((uint64)current_ptr_5 == (uint64)(top_of_action_table + action_matrix_size_5));
    return result;
}

void FiniteDimensionalModule_setActionVector(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_idx,
    int input_degree, uint input_idx,
    Vector *output
){
    input_degree -= module->module.min_degree;
    uint output_degree = input_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector *output_vector = module->actions[input_degree][output_degree][operation_idx][input_idx];
    Vector_assign(output_vector, output);
}

void FiniteDimensionalModule_setAction(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_idx,
    int input_degree, uint input_idx,
    uint *output
){
    // printf("    operation_degree: %d, operation_idx: %d, input_degree: %d, input_idx: %d\n", operation_degree, operation_idx, input_degree, input_idx);
    input_degree -= module->module.min_degree;
    assert(input_degree >= 0);
    uint output_degree = input_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector *output_vector = module->actions[input_degree][output_degree][operation_idx][input_idx];  
    // array_print("    output: %s\n", output, output_vector->dimension);
    Vector_pack(output_vector, output);
}


Vector *FiniteDimensionalModule_getAction(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_index,
    int module_degree, uint module_index
){
    module_degree -= module->module.min_degree;
    int output_degree = module_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector *output_vector = module->actions[module_degree][output_degree][operation_index][module_index];
    return output_vector;
}

bool FiniteDimensionalModule_computeBasis(Module *this __attribute__((unused)), int degree __attribute__((unused))){
    return true;
}

uint FiniteDimensionalModule_getDimension(Module *this, int degree){
    assert(degree >= this->min_degree);
    FiniteDimensionalModule *module = (FiniteDimensionalModule *) this;
    if(degree >= module->max_basis_degree ){
        return 0;
    }
    return module->graded_dimension[degree - this->min_degree];
    return 0;
}

void FiniteDimensionalModule_actOnBasis(Module *this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index){
    FiniteDimensionalModule *module = ((FiniteDimensionalModule*)this);
    assert(op_index < algebra_getDimension(this->algebra, op_degree, mod_degree));
    assert(mod_index < module_getDimension(this, mod_degree));
    uint output_dimension = module_getDimension(this, mod_degree + op_degree);    
    if(mod_degree + op_degree >= module->max_basis_degree || output_dimension == 0){
        return;
    }  
    // Why do we take this slice?
    char output_block_memory[Vector_getContainerSize(this->p)];    
    Vector *output_block = Vector_initialize(this->p, output_block_memory, NULL, 0, 0);     
    Vector_slice(output_block, result, 0, output_dimension); 
    Vector *output = FiniteDimensionalModule_getAction(module, op_degree, op_index, mod_degree, mod_index);
    Vector_add(output_block, output, coeff);
}

typedef struct {
    Module module;
    uint *number_of_generators_in_degree;
    // private fields
    uint computed_degree;
    FreeModuleOperationGeneratorPair **basis_element_to_opgen_table;
    uint ***generator_to_index_table;
} FreeModuleInternal;

FreeModule *FreeModule_construct(Algebra *algebra, int min_degree, int max_degree){
    uint num_degrees = max_degree - min_degree;
    size_t module_size =  sizeof(FreeModuleInternal) 
        + num_degrees * sizeof(uint) // number_of_generators_in_degree
        + num_degrees * sizeof(FreeModuleOperationGeneratorPair*) //basis_element_to_opgen_table
        + num_degrees * sizeof(uint **) // generator_to_index_table
    ;
    FreeModuleInternal *module = malloc(module_size);
    memset(module, 0, module_size);
    module->module.algebra = algebra;
    module->module.p = algebra->p;
    module->module.computeBasis = FreeModule_computeBasis;
    module->module.getDimension = FreeModule_getDimension;
    module->module.actOnBasis = FreeModule_actOnBasis;
    module->module.max_degree = max_degree;
    module->module.min_degree = min_degree;
    module->computed_degree = 0;
    module->number_of_generators_in_degree = (uint *)(module + 1);
    module->basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair**)(
        module->number_of_generators_in_degree + num_degrees
    );    
    module->generator_to_index_table = (uint ***)(module->basis_element_to_opgen_table + num_degrees);
    assert((char*)(module->generator_to_index_table + num_degrees) == (char*)module + module_size);
    return (FreeModule*)module;
}

void FreeModule_free(FreeModule *module){
    if(module == NULL){
        return;
    }
    FreeModuleInternal *M = (FreeModuleInternal*) module;
    for(uint i = 0; i < module->module.max_degree - module->module.min_degree; i++){
        free(M->generator_to_index_table[i]);
    }
    free(module);
}

bool FreeModule_computeBasis(Module *this __attribute__((unused)), int degree __attribute__((unused))){
    return true;
}

uint FreeModule_getDimension(Module *this, int degree){
    assert(degree >= this->min_degree);
    assert(degree < this->max_degree);
    FreeModule *module = (FreeModule*) this;
    uint result = 0;
    for(int i = this->min_degree; i <= degree; i++){
        // for the excess make sure to use math degree
        result += FreeModule_getNumberOfGensInDegree(module, i)
            * algebra_getDimension(this->algebra, degree - i, i);
    }
    return result;
}

// Run FreeModule_ConstructBlockOffsetTable(module, degree) before using this on an input in that degree
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, int op_deg, uint op_idx, int module_degree, uint module_idx){
    assert(op_idx < algebra_getDimension(this->algebra, op_deg, module_degree));
    assert(FreeModule_getDimension(this, op_deg + module_degree) <= result->dimension);
    FreeModuleInternal *module = (FreeModuleInternal *) this;
    FreeModuleOperationGeneratorPair operation_generator = FreeModule_indexToOpGen((FreeModule*)module, module_degree, module_idx);
    int module_operation_degree = operation_generator.operation_degree;
    uint module_operation_index = operation_generator.operation_index;
    int generator_degree = operation_generator.generator_degree; 
    uint generator_index  = operation_generator.generator_index;
    // Now all of the output elements are going to be of the form s * x. Find where such things go in the output vector.
    uint num_ops = algebra_getDimension(this->algebra, module_operation_degree + op_deg, generator_degree);
    uint output_block_min = FreeModule_operationGeneratorToIndex((FreeModule*)module, module_operation_degree + op_deg, 0, generator_degree, generator_index);

    uint output_block_max = output_block_min + num_ops;
    char output_block_memory[Vector_getContainerSize(this->p)];    
    Vector *output_block = Vector_initialize(this->p, output_block_memory, NULL, 0, 0);     
    Vector_slice(output_block, result, output_block_min, output_block_max); 
    // Now we multiply s * r and write the result to the appropriate position.
    algebra_multiplyBasisElements(module->module.algebra, output_block, coeff, op_deg, op_idx, module_operation_degree, module_operation_index, generator_degree);
}


// Compute tables:
//    basis element index     --> operator, generator pair
//    a generator  --> where does that generator's block start?
void FreeModule_ConstructBlockOffsetTable(FreeModule *M, int degree){
    assert(degree>=M->module.min_degree);
    assert(degree < M->module.max_degree);
    int shifted_degree = degree - M->module.min_degree;
    FreeModuleInternal *module = (FreeModuleInternal *) M;
    if(module->generator_to_index_table[shifted_degree] != NULL){
        return;
    }
    // gen_to_idx goes gen_degree => gen_idx => start of block.
    // so gen_to_idx_size should be (number of possible degrees + 1) * sizeof(uint*) + number of gens * sizeof(uint).
    // The other part of the table goes idx => opgen
    // The size should be (number of basis elements in current degree) * sizeof(FreeModuleOperationGeneratorPair)
    // A basis element in degree n comes from a generator in degree i paired with an operation in degree n - i.
    size_t gen_to_idx_size = (degree - M->module.min_degree + 1) * sizeof(uint*);
    size_t total_size = 0;
    for(int gen_deg = M->module.min_degree; gen_deg <= degree; gen_deg++){
        int op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(module->module.algebra, op_deg, gen_deg);
        uint num_gens = FreeModule_getNumberOfGensInDegree(M, gen_deg);
        gen_to_idx_size += num_gens * sizeof(uint);
        total_size += num_gens * num_ops * sizeof(FreeModuleOperationGeneratorPair);
    }
    total_size += gen_to_idx_size;
    char *memory = realloc(module->generator_to_index_table[shifted_degree], total_size);
    uint **generator_to_index_table = (uint**) memory;
    FreeModuleOperationGeneratorPair * basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair *)(memory + gen_to_idx_size);
    module->generator_to_index_table[shifted_degree] = generator_to_index_table;
    module->basis_element_to_opgen_table[shifted_degree] = basis_element_to_opgen_table;
    uint **gentoidx_degree_ptr = generator_to_index_table;
    uint *gentoidx_index_ptr = (uint*)(generator_to_index_table + (shifted_degree + 1));
    uint offset = 0;
    uint generator = 0;
    for(int gen_deg = M->module.min_degree; gen_deg <= degree; gen_deg++){
        *gentoidx_degree_ptr = gentoidx_index_ptr;
        int op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(module->module.algebra, op_deg, gen_deg);
        for(uint gen_idx = 0; gen_idx < FreeModule_getNumberOfGensInDegree(M, gen_deg); gen_idx++){
            *gentoidx_index_ptr = offset;
            for(uint op_idx = 0; op_idx < num_ops; op_idx++){
                basis_element_to_opgen_table->generator_degree = gen_deg;
                basis_element_to_opgen_table->generator_index  = gen_idx;
                basis_element_to_opgen_table->operation_degree = op_deg;
                basis_element_to_opgen_table->operation_index  = op_idx;
                basis_element_to_opgen_table ++;
            }
            gentoidx_index_ptr++;
            generator ++;
            offset += num_ops;
        }
        gentoidx_degree_ptr++;
    }
    assert(gentoidx_degree_ptr == generator_to_index_table + (shifted_degree + 1));
    assert((char*)gentoidx_index_ptr == memory + gen_to_idx_size);
    assert((char*)basis_element_to_opgen_table == memory + total_size);
}



uint FreeModule_getNumberOfGensInDegree(FreeModule *this, int degree){
    return this->number_of_generators_in_degree[degree - this->module.min_degree];
}

uint FreeModule_operationGeneratorToIndex(FreeModule *this, int op_deg, uint op_idx, int gen_deg, uint gen_idx){
    FreeModuleInternal *module = (FreeModuleInternal *)this;
    gen_deg -= this->module.min_degree;
    assert(gen_deg >= 0);
    int deg = op_deg + gen_deg;
    assert(module->generator_to_index_table[deg]!=NULL);
    uint block_idx = module->generator_to_index_table[deg][gen_deg][gen_idx];
    return block_idx + op_idx;
}

FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, int degree, uint index){
    degree -= this->module.min_degree;
    FreeModuleInternal *module = (FreeModuleInternal*) this;
    return module->basis_element_to_opgen_table[degree][index];
}



// FreeModuleHomomorphisms

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
    FreeModuleInternal *module = (FreeModuleInternal*) f->source;
    uint p = module->module.p;
    uint dimension = module_getDimension(f->target, degree);
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

Vector *FreeModuleHomomorphism_getOutput(FreeModuleHomomorphism *f, int input_degree, uint input_index){
    assert(input_degree - f->source->module.min_degree >= 0);
    assert(input_degree < f->max_computed_degree);
    assert(input_index < module_getDimension((Module*)f->source, input_degree));
    return f->outputs[input_degree - f->source->module.min_degree][input_index];
}

void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, int input_degree, uint input_index, Vector *output){
    assert(output->dimension == module_getDimension(f->target, input_degree));
    assert(output->offset == 0);
    assert(input_index < FreeModule_getNumberOfGensInDegree(f->source, input_degree));
    Vector_assign(f->outputs[input_degree - f->source->module.min_degree][input_index], output);
}

// Run FreeModule_ConstructBlockOffsetTable(source, degree) before using this on an input in that degree
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, int input_degree, uint input_index){
    int shifted_input_degree = input_degree - f->source->module.min_degree;
    assert(shifted_input_degree > 0);
    assert(input_degree < f->max_degree);
    assert(((FreeModuleInternal*)f->source)->basis_element_to_opgen_table[shifted_input_degree] != NULL);
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
            module_actOnBasis(f->target, result, c, operation_degree, operation_index, generator_degree, it.index);
        }
    }
}


// result should be big enough to hold output (how big is that?)
// Well it should have dim(target) columns and dim(source) rows. I guess this is 
// the transpose of the usual convention.
void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, int degree){
    assert(degree < f->max_degree);
    assert(module_getDimension(&f->source->module, degree) <= result->rows);
    assert(module_getDimension(f->target, degree) <= result->columns);
    // The shorter implementation if we do FreeModuleConstructBlockOffsetTable first.
    // Maybe we ought to do that...
    // for(int i = 0; i < module_getDimension(&f->source->module, degree); i++){
    //     FreeModuleHomomorphism_applyToBasisElement(f, result->matrix[i], 1, degree, i);
    // }    
    FreeModule *source = f->source;
    Algebra *algebra = source->module.algebra;
    // 
    uint i = 0;
    for(int gen_deg = f->source->module.min_degree; gen_deg <= degree; gen_deg++){
        int op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(algebra, op_deg, gen_deg);
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
                        module_actOnBasis(f->target, result->matrix[i], it.value, op_deg, op_idx, gen_deg, it.index);
                    }
                }
                i++;                
            }
        }
    }
}

Kernel *Kernel_construct(uint p, uint rows, uint columns){
    assert(columns < MAX_DIMENSION);
    Kernel *k = malloc(
        sizeof(Kernel) 
        + columns * sizeof(uint)
        + Matrix_getSize(p, rows, columns) * sizeof(uint64)
    );
    k->column_to_pivot_row = (int*)(k + 1);
    k->kernel = Matrix_initialize((char*)(k->column_to_pivot_row + columns), p, rows, columns);
    return k;
}

void Kernel_free(Kernel *k){
    free(k);
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