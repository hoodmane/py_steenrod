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
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint max_generator_degree, uint *graded_dimension);

FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, uint max_degree, uint *graded_dimension){
    max_degree ++;
    array_print(graded_dimension, max_degree);
    FiniteDimensionalModule *result = FiniteDimensionalModule_allocate(algebra, max_degree, graded_dimension);
    result->module.p = algebra->p;
    result->module.algebra = algebra;
    result->module.computeBasis = FiniteDimensionalModule_computeBasis;
    result->module.getDimension = FiniteDimensionalModule_getDimension;
    result->module.actOnBasis = FiniteDimensionalModule_actOnBasis;
    result->max_degree = max_degree;
    result->dimension = 0;
    for(uint i = 0; i <= max_degree; i++){
        result->dimension += graded_dimension[i];
    }
    memcpy(result->graded_dimension, graded_dimension, max_degree * sizeof(uint));
    return result;
}

void FiniteDimensionalModule_free(FiniteDimensionalModule *module){
    free(module);
}

// This is the grossest allocator.
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint max_generator_degree, uint *graded_dimension){
    uint p = algebra->p;
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
    for(uint degree = 0; degree < max_generator_degree; degree++) {
        if(graded_dimension[degree] != 0){
            number_of_nonempty_degrees ++;
        }
    }
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> (out_index) -> value
    //  ****    -> ***       -> **Vector   -> *Vector    -> Vector -> uint
    action_matrix_size_1 += max_generator_degree * sizeof(Vector***);
    action_matrix_size_2 += number_of_nonempty_degrees * max_generator_degree * sizeof(Vector**);
    for(uint input_degree = 0; input_degree < max_generator_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            continue;
        }
        for(uint output_degree = input_degree + 1; output_degree < max_generator_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                continue;
            }
            uint number_of_operations = algebra_getDimension(algebra, output_degree - input_degree);
            action_matrix_size_3 += sizeof(Vector*) * number_of_operations;
            action_matrix_size_4 += Vector_getContainerSize(p) * number_of_operations * graded_dimension[input_degree];
            uint vectorSize = Vector_getSize(p, graded_dimension[output_degree], 0);
            action_matrix_size_5 += number_of_operations * graded_dimension[input_degree] * vectorSize;
        }
    }

    action_matrix_size_2 += action_matrix_size_1;
    action_matrix_size_3 += action_matrix_size_2;
    action_matrix_size_4 += action_matrix_size_3;
    action_matrix_size_5 += action_matrix_size_4;

    FiniteDimensionalModule *result = malloc(
            sizeof(FiniteDimensionalModule)
            + max_generator_degree * sizeof(uint)
            + action_matrix_size_5
    );
    result->graded_dimension = (uint *) (result + 1);
    char *top_of_action_table = (char *) (result + 1) + max_generator_degree *sizeof(uint);
    Vector ****current_ptr_1 = (Vector ****) top_of_action_table;
    Vector ***current_ptr_2 = (Vector ***) (top_of_action_table + action_matrix_size_1);
    Vector **current_ptr_3 = (Vector **) (top_of_action_table + action_matrix_size_2);
    Vector *current_ptr_4 = (Vector *) (top_of_action_table + action_matrix_size_3);
    char *current_ptr_5 = (char *) (top_of_action_table + action_matrix_size_4);
    result->actions = current_ptr_1; 
    memset(top_of_action_table, 0, action_matrix_size_5);
    for(int input_degree = 0; input_degree < max_generator_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            current_ptr_1 ++;
            continue;
        }
        *current_ptr_1 = current_ptr_2;
        current_ptr_2 += input_degree + 1;
        for(int output_degree = input_degree + 1; output_degree < max_generator_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                current_ptr_2 ++;
                continue;
            }
            *current_ptr_2 = current_ptr_3;
            uint vectorSize = Vector_getSize(p, graded_dimension[output_degree], 0);
            uint number_of_operations = algebra_getDimension(algebra, output_degree - input_degree);
            for(int operation_idx = 0; operation_idx < number_of_operations; operation_idx ++){
                *current_ptr_3 = current_ptr_4;
                for(int input_idx = 0; input_idx < graded_dimension[input_degree]; input_idx ++ ){
                    Vector_initialize(p, (char*)current_ptr_4, current_ptr_5, graded_dimension[output_degree], 0);
                    // ... gross:
                    current_ptr_4 = (Vector*)(((char*)current_ptr_4) + Vector_getContainerSize(p));
                    current_ptr_5 += graded_dimension[output_degree] * vectorSize;
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
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    Vector *output
){
    uint output_degree = input_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector *output_vector = &module->actions[input_degree][output_degree][operation_idx][input_idx];
    Vector_assign(output_vector, output);
}

void FiniteDimensionalModule_setAction(
    FiniteDimensionalModule *module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    uint *output
){
    uint output_degree = input_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector *output_vector = &module->actions[input_degree][output_degree][operation_idx][input_idx];
    Vector_pack(output_vector, output);
}

bool FiniteDimensionalModule_computeBasis(Module *this, uint dimension){
    return true;
}

uint FiniteDimensionalModule_getDimension(Module *this, uint degree){
    FiniteDimensionalModule *module = (FiniteDimensionalModule *) this;
    if(degree < module->max_degree ){
        return module->graded_dimension[degree];
    }
    return 0;
}

void FiniteDimensionalModule_actOnBasis(Module *this, Vector *result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_index){
    FiniteDimensionalModule *module = ((FiniteDimensionalModule*)this);
    assert(op_index < algebra_getDimension(this->algebra, op_degree));
    assert(mod_index < module->graded_dimension[mod_degree]);    
    if(mod_degree + op_degree >= module->max_degree || module->graded_dimension[mod_degree + op_degree] == 0){
        return;
    }
    char output_block_memory[Vector_getContainerSize(this->p)];    
    Vector *output_block = (Vector*)output_block_memory;     
    Vector_slice(output_block, result, 0, module->graded_dimension[mod_degree + op_degree]); 
    Vector_add(output_block, &module->actions[mod_degree][mod_degree + op_degree][op_index][mod_index], coeff);
}

typedef struct {
    Module module;
    uint max_generator_degree;
    uint max_degree;
    uint number_of_generators;
    uint *number_of_generators_in_degree;
    // private fields
    uint computed_degree;
    FreeModuleOperationGeneratorPair **basis_element_to_opgen_table;
    uint ***generator_to_index_table;
} FreeModuleInternal;

FreeModule *FreeModule_construct(Algebra *algebra, uint max_generator_degree, uint max_degree){
    FreeModuleInternal *module = malloc(
        sizeof(FreeModuleInternal) 
        + (max_degree + 1) * sizeof(uint) // number_of_generators_in_degree
        + (max_degree + 1) * sizeof(FreeModuleOperationGeneratorPair*) //basis_element_to_opgen_table
        + (max_degree + 1) * sizeof(uint **) // generator_to_index_table
    );
    memset(module + 1, 0, 
        (max_degree + 1) * sizeof(uint) // number_of_generators_in_degree
        + (max_degree + 1) * sizeof(FreeModuleOperationGeneratorPair*) //basis_element_to_opgen_table
        + (max_degree + 1) * sizeof(uint **) // generator_to_index_table
    );
    module->module.algebra = algebra;
    module->module.p = algebra->p;
    module->module.computeBasis = FreeModule_computeBasis;
    module->module.getDimension = FreeModule_getDimension;
    module->module.actOnBasis = FreeModule_actOnBasis;
    module->max_generator_degree = max_generator_degree;
    module->max_degree = max_degree;
    module->computed_degree = 0;
    module->number_of_generators = 0;
    module->number_of_generators_in_degree = (uint *)(module + 1);
    module->generator_to_index_table = (uint ***)(module->number_of_generators_in_degree + (max_degree + 1));
    module->basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair**)(
        module->generator_to_index_table + (max_degree + 1)
    );
    memset(module->number_of_generators_in_degree, 0, (max_generator_degree + 1) * sizeof(uint));
    return (FreeModule*)module;
}

void FreeModule_free(FreeModule *module){
    // TODO: Free other parts of the structure.
    free(module);
}

bool FreeModule_computeBasis(Module *this, uint degree){
    return true;
}

uint FreeModule_getDimension(Module *this, uint degree){
    FreeModule *module = (FreeModule*) this;
    uint result = 0;
    for(int i = 0; i <= module->max_generator_degree && i<= degree; i++){
        result += module->number_of_generators_in_degree[i] * algebra_getDimension(this->algebra, degree - i);
    }
    return result;
}

// Run FreeModule_ConstructBlockOffsetTable(module, degree) before using this on an input in that degree
void FreeModule_actOnBasis(Module *this, Vector *result, uint coeff, uint op_deg, uint op_idx, uint module_degree, uint module_idx){
    assert(op_idx < algebra_getDimension(this->algebra, op_deg));
    assert(FreeModule_getDimension(this, op_deg + module_degree) <= result->dimension);
    FreeModuleInternal *module = (FreeModuleInternal *) this;
    FreeModuleOperationGeneratorPair operation_generator = module->basis_element_to_opgen_table[module_degree][module_idx];
    uint module_operation_degree = operation_generator.operation_degree;
    uint module_operation_index = operation_generator.operation_index;
    uint generator_degree = operation_generator.generator_degree; 
    uint generator_index  = operation_generator.generator_index;
    // Now all of the output elements are going to be of the form s * x. Find where such things go in the output vector.
    uint num_ops = algebra_getDimension(this->algebra, module_operation_degree + op_deg);
    uint output_block_min = module->generator_to_index_table[module_degree + op_deg][generator_degree][generator_index];
    uint output_block_max = output_block_min + num_ops;
    char output_block_memory[Vector_getContainerSize(this->p)];    
    Vector *output_block = (Vector*)output_block_memory;     
    Vector_slice(output_block, result, output_block_min, output_block_max); 
    // Now we multiply s * r and write the result to the appropriate position.
    algebra_multiplyBasisElements(module->module.algebra, output_block, coeff, op_deg, op_idx, module_operation_degree, module_operation_index);
}

// Compute tables:
//    basis element index     --> operator, generator pair
//    a generator  --> where does that generator's block start?
void FreeModule_ConstructBlockOffsetTable(FreeModule *M, uint degree){
    FreeModuleInternal *module = (FreeModuleInternal *) M;
    if(module->generator_to_index_table[degree] != NULL){
        return;
    }
    uint max_gen_degree = module->max_generator_degree < degree ? module->max_generator_degree : degree;
    size_t gen_to_idx_size = (max_gen_degree + 1) * sizeof(uint*);
    size_t total_size = 0;
    for(uint gen_deg = 0; gen_deg <= max_gen_degree; gen_deg++){
        uint op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(module->module.algebra, op_deg);
        uint num_gens = module->number_of_generators_in_degree[gen_deg];
        gen_to_idx_size += num_gens * sizeof(uint);
        total_size += num_gens * num_ops * sizeof(FreeModuleOperationGeneratorPair);
    }
    total_size += gen_to_idx_size;
    char *memory = realloc(module->generator_to_index_table[degree], total_size);
    uint **generator_to_index_table = (uint**) memory;
    FreeModuleOperationGeneratorPair * basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair *)(memory + gen_to_idx_size);
    module->generator_to_index_table[degree] = generator_to_index_table;
    module->basis_element_to_opgen_table[degree] = basis_element_to_opgen_table;
    uint **gentoidx_degree_ptr = generator_to_index_table;
    uint *gentoidx_index_ptr = (uint*)(generator_to_index_table + max_gen_degree + 1);
    uint offset = 0;
    uint generator = 0;
    for(uint gen_deg = 0; gen_deg <= max_gen_degree; gen_deg++){
        *gentoidx_degree_ptr = gentoidx_index_ptr;
        uint op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(module->module.algebra, op_deg);
        for(uint gen_idx = 0; gen_idx < module->number_of_generators_in_degree[gen_deg]; gen_idx++){
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
    assert(gentoidx_degree_ptr == generator_to_index_table + max_gen_degree + 1);
    assert((char*)gentoidx_index_ptr == memory + gen_to_idx_size);
    assert((char*)basis_element_to_opgen_table == memory + total_size);
}

uint FreeModule_operationGeneratorToIndex(FreeModule *this, uint op_deg, uint op_idx, uint gen_deg, uint gen_idx){
    FreeModuleInternal *module = (FreeModuleInternal *)this;
    uint deg = op_deg + gen_deg;
    assert(module->generator_to_index_table[deg]!=NULL);
    uint block_idx = module->generator_to_index_table[deg][gen_deg][gen_idx];
    return block_idx + op_idx;
}

FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, uint degree, uint index){
    FreeModuleInternal *module = (FreeModuleInternal*) this;
    return module->basis_element_to_opgen_table[degree][index];
}



// FreeModuleHomomorphisms

FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule *source, Module *target, uint max_degree){
    FreeModuleHomomorphism *f = malloc(
        sizeof(FreeModuleHomomorphism) 
        + (max_degree + 1) * sizeof(Vector**) // outputs
        + (max_degree + 1) * sizeof(Matrix*)  // coimage_to_image_iso
        + (max_degree + 1) * sizeof(Kernel*)  // kernel
    );
    f->source = source;
    f->target = target;
    f->max_computed_degree = -1;
    f->outputs = (Vector***)(f + 1);
    f->coimage_to_image_isomorphism = (Matrix**)(f->outputs + max_degree + 1);
    f->kernel = (Kernel**)(f->coimage_to_image_isomorphism + max_degree + 1);
    memset(f+1, 0,(max_degree + 1)*(sizeof(Vector**) + sizeof(Matrix*) + sizeof(Kernel*)));
    return f;
}

void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism *f, uint degree, uint num_gens){
    FreeModuleInternal *module = (FreeModuleInternal*) f->source;
    uint p = module->module.p;
    uint dimension = module_getDimension(f->target, degree);
    uint vector_size = Vector_getSize(p, dimension, 0);
    f->outputs[degree] = (Vector**)malloc(
        num_gens * sizeof(Vector*) 
        + num_gens * Vector_getContainerSize(p)
        + num_gens * vector_size
    );
    Vector **vector_ptr_ptr = f->outputs[degree];
    char *vector_container_ptr = (char*)(vector_ptr_ptr + num_gens);
    char *vector_memory_ptr = vector_container_ptr + num_gens * Vector_getContainerSize(p);
    for(uint i = 0; i < num_gens; i++){
        f->outputs[degree][i] = Vector_initialize(p, vector_container_ptr, vector_memory_ptr, dimension, 0);
        vector_ptr_ptr ++;
        vector_container_ptr += Vector_getContainerSize(p);
        vector_memory_ptr += vector_size;
    }
    assert(vector_ptr_ptr == f->outputs[degree] + num_gens);
    assert(vector_container_ptr == (char*)(vector_ptr_ptr) + num_gens * Vector_getContainerSize(p));
    assert(vector_memory_ptr == vector_container_ptr + num_gens * vector_size);

}

void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, uint input_degree, uint input_index, Vector *output){
    assert(output->dimension == module_getDimension(f->target, input_degree));
    assert(output->offset == 0);
    assert(input_index < f->source->number_of_generators_in_degree[input_degree]);
    Vector_assign(f->outputs[input_degree][input_index], output);
}

// Run FreeModule_ConstructBlockOffsetTable(source, degree) before using this on an input in that degree
void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism *f, Vector *result, uint coeff, uint input_degree, uint input_index){
    assert(((FreeModuleInternal*)f->source)->basis_element_to_opgen_table[input_degree] != NULL);
    assert(input_index < FreeModule_getDimension((Module*)f->source, input_degree));
    FreeModuleOperationGeneratorPair operation_generator = 
        ((FreeModuleInternal*)f->source)->basis_element_to_opgen_table[input_degree][input_index];
    uint operation_degree = operation_generator.operation_degree;
    uint operation_index = operation_generator.operation_index;
    uint generator_degree = operation_generator.generator_degree;
    uint generator_index = operation_generator.generator_index;
    Vector *output_on_generator = f->outputs[generator_degree][generator_index];
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
void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism *f, Matrix *result, uint degree){
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
    uint max_degree = source->max_generator_degree < degree ? source->max_generator_degree : degree;
    for(uint gen_deg = 0; gen_deg <= max_degree; gen_deg++){
        uint op_deg = degree - gen_deg;
        uint num_ops = algebra_getDimension(algebra, op_deg);
        for(uint gen_idx = 0; gen_idx < source->number_of_generators_in_degree[gen_deg]; gen_idx++){
            for(uint op_idx = 0; op_idx < num_ops; op_idx++){
                Vector *output_on_generator = f->outputs[gen_deg][gen_idx];
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