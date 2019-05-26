//
// Created by Hood on 5/20/2019.
//

#include <stdlib.h>
#include <assert.h>

#include "modules.h"
#include "FpVector.h"
#include "algebra.h"

bool FiniteDimensionalModule_compute_basis(Module *this, uint dimension);
uint FiniteDimensionalModule_get_dimension(Module* this, uint degree);
void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_index);

// The allocator is horrendous so we're going to separate it out.
FiniteDimensionalModule * allocateFiniteDimensionalModule(Algebra * algebra, uint max_generator_degree, uint * number_of_generators_in_degree);

FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, uint max_generator_degree, uint * number_of_generators_in_degree){
    max_generator_degree ++;
    FiniteDimensionalModule * result = allocateFiniteDimensionalModule(algebra, max_generator_degree, number_of_generators_in_degree);
    result->module.p = algebra->p;
    result->module.algebra = algebra;
    result->module.compute_basis = FiniteDimensionalModule_compute_basis;
    result->module.get_dimension = FiniteDimensionalModule_get_dimension;
    result->module.act_on_basis = FiniteDimensionalModule_act_on_basis;
    result->max_degree = max_generator_degree;
    result->dimension = 0;
    for(uint i = 0; i < max_generator_degree; i++){
        result->dimension += number_of_generators_in_degree[i];
    }
    memcpy(result->number_of_basis_elements_in_degree, number_of_generators_in_degree, max_generator_degree * sizeof(uint));
    return result;
}

void freeFiniteDimensionalModule(FiniteDimensionalModule * module){
    free(module);
}

// This is the grossest allocator.
FiniteDimensionalModule * allocateFiniteDimensionalModule(Algebra * algebra, uint max_generator_degree, uint * number_of_generators_in_degree){
    uint p = algebra->p;
    VectorInterface vectorInterface = algebra->vectorInterface;
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
        if(number_of_generators_in_degree[degree] != 0){
            number_of_nonempty_degrees ++;
        }
    }
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> (out_index) -> value
    //  ****    -> ***       -> **Vector   -> *Vector    -> Vector -> uint
    action_matrix_size_1 += max_generator_degree * sizeof(Vector***);
    action_matrix_size_2 += number_of_nonempty_degrees * max_generator_degree * sizeof(Vector**);
    for(int input_degree = 0; input_degree < max_generator_degree; input_degree++){
        if(number_of_generators_in_degree[input_degree] == 0){
            continue;
        }
        for(int output_degree = input_degree + 1; output_degree < max_generator_degree; output_degree++){
            if(number_of_generators_in_degree[output_degree] == 0){
                continue;
            }
            uint number_of_operations = algebra_get_dimension(algebra, output_degree - input_degree);
            action_matrix_size_3 += sizeof(Vector*)       * number_of_operations;
            action_matrix_size_4 += sizeof(Vector)        * number_of_operations * number_of_generators_in_degree[input_degree];
            uint vectorSize = vectorInterface.getSize(p, number_of_generators_in_degree[output_degree], 0);
            action_matrix_size_5 += sizeof(uint64) * number_of_operations * number_of_generators_in_degree[input_degree] * vectorSize;
        }
    }

    action_matrix_size_2 += action_matrix_size_1;
    action_matrix_size_3 += action_matrix_size_2;
    action_matrix_size_4 += action_matrix_size_3;
    action_matrix_size_5 += action_matrix_size_4;

    FiniteDimensionalModule * result = malloc(
            sizeof(FiniteDimensionalModule)
            + max_generator_degree * sizeof(uint)
            + action_matrix_size_5
    );
    result->number_of_basis_elements_in_degree = (uint * ) (result + 1);
    void * top_of_action_table = (void *) (result + 1) + max_generator_degree * sizeof(uint);
    Vector **** current_ptr_1 = (Vector ****) top_of_action_table;
    Vector *** current_ptr_2 = (Vector ***) (top_of_action_table + action_matrix_size_1);
    Vector ** current_ptr_3 = (Vector **) (top_of_action_table + action_matrix_size_2);
    Vector * current_ptr_4 = (Vector *) (top_of_action_table + action_matrix_size_3);
    uint64 * current_ptr_5 = (uint64 *) (top_of_action_table + action_matrix_size_4);
    result->actions = top_of_action_table;
    memset(top_of_action_table, 0, action_matrix_size_5);
    for(int input_degree = 0; input_degree < max_generator_degree; input_degree++){
        if(number_of_generators_in_degree[input_degree] == 0){
            continue;
        }
        *current_ptr_1 = current_ptr_2;
        current_ptr_2 += input_degree + 1;
        for(int output_degree = input_degree + 1; output_degree < max_generator_degree; output_degree++){
            if(number_of_generators_in_degree[output_degree] == 0){
                continue;
            }
            *current_ptr_2 = current_ptr_3;
            uint number_of_operations = algebra_get_dimension(algebra, output_degree - input_degree);
            for(int operation_idx = 0; operation_idx < number_of_operations; operation_idx ++){
                *current_ptr_3 = current_ptr_4;
                for(int input_idx = 0; input_idx < number_of_generators_in_degree[input_degree]; input_idx ++ ){
                    vectorInterface.initialize(p, (uint64*)current_ptr_4, current_ptr_5, number_of_generators_in_degree[output_degree], 0);
                    current_ptr_4 = (Vector*)(((uint64*)current_ptr_4) + vectorInterface.container_size);
                    current_ptr_5 += number_of_generators_in_degree[output_degree];
                }
                current_ptr_3 ++;
            }
            current_ptr_2 ++;
        }
        current_ptr_1 ++;
    }
//    assert((uint)current_ptr_1 == (uint)(top_of_action_table + action_matrix_size_1));
//    assert((uint)current_ptr_2 == (uint)(top_of_action_table + action_matrix_size_2));
//    assert((uint)current_ptr_3 == (uint)(top_of_action_table + action_matrix_size_3));
//    assert((uint)current_ptr_4 == (uint)(top_of_action_table + action_matrix_size_4));
//    assert((uint)current_ptr_5 == (uint)(top_of_action_table + action_matrix_size_5));
    return result;
}

void addActionToFiniteDimensionalModule(
    FiniteDimensionalModule * module,
    uint operation_degree, uint operation_idx,
    uint input_degree, uint input_idx,
    uint * output
){
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector * output_vector = &module->actions[input_degree][input_degree + operation_degree][operation_idx][input_idx];
    VectorInterface vectorInterface = module->module.algebra->vectorInterface;
    vectorInterface.pack(output_vector, output);
}

bool FiniteDimensionalModule_compute_basis(Module *this, uint dimension){
    return true;
}

uint FiniteDimensionalModule_get_dimension(Module* this, uint degree){
    FiniteDimensionalModule * module = (FiniteDimensionalModule *) this;
    if(degree <= module->max_degree ){
        return module->number_of_basis_elements_in_degree[degree];
    }
    return 0;
}

void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_index){
    FiniteDimensionalModule* module = ((FiniteDimensionalModule*)this);
    VectorInterface vectorInterface = this->algebra->vectorInterface;
    vectorInterface.add(result, &module->actions[mod_degree][mod_degree + op_degree][op_index][mod_index], coeff);
}


bool FreeModule_compute_basis(Module* this, uint degree);
uint FreeModule_get_dimension(Module* this, uint degree);
void FreeModule_act_on_basis(Module * this, Vector * result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_idx);

FreeModule * constructFreeModule(Algebra * algebra, uint max_degree){
    FreeModule * module = malloc(sizeof(FreeModule) + 2*(max_degree + 1) * sizeof(uint **));
    module->module.algebra = algebra;
    module->module.p = algebra->p;
    module->module.compute_basis = FreeModule_compute_basis;
    module->module.get_dimension = FreeModule_get_dimension;
    module->module.act_on_basis = FreeModule_act_on_basis;
    module->max_degree = max_degree;
    module->computed_degree++;
    module->generator_to_index_table = (uint **)(module + 1);
    module->basis_element_to_opgen_table = (FreeModuleOperationGeneratorPair**)(module->generator_to_index_table + (max_degree + 1));
    return module;
}

void freeFreeModule(FreeModule * module){
    free(module);
}

bool FreeModule_compute_basis(Module* this, uint degree){
    return true;
}


// Compute tables:
//    basis element index     --> operator, generator pair
//    a generator  --> where does that generator's block start?
void FreeModuleComputeBlockOffsetTable(FreeModule * module, uint * generator_to_index_table, FreeModuleOperationGeneratorPair * basis_element_to_opgen_table, uint degree){
    uint offset = 0;
    uint generator = 0;
    for(uint gen_deg = 0; gen_deg <= module->max_generator_degree && gen_deg < degree; gen_deg++){
        uint op_deg = degree - gen_deg;
        uint num_ops = algebra_get_dimension(module->module.algebra, op_deg);
        for(uint gen_idx = 0; gen_idx < module->number_of_generators_in_degree[gen_deg]; gen_idx++){
            *generator_to_index_table = offset;
            for(uint op_idx = 0; op_idx < num_ops; op_idx++){
                basis_element_to_opgen_table->generator_degree = gen_deg;
                basis_element_to_opgen_table->generator_index  = gen_idx;
                basis_element_to_opgen_table->operation_degree = op_deg;
                basis_element_to_opgen_table->operation_index  = op_idx;
                basis_element_to_opgen_table ++;
            }
            generator_to_index_table ++;
            generator ++;
            offset += num_ops;
        }
    }
}

uint FreeModule_get_dimension(Module* this, uint degree){
    FreeModule * module = (FreeModule*) this;
    uint result = 0;
    for(int i = 0; i <= degree; i++){
        result += module->number_of_generators_in_degree[i] * algebra_get_dimension(this->algebra, degree - i);
    }
    return result;
}


void FreeModule_act_on_basis(Module * this, Vector * result, uint coeff, uint op_deg, uint op_idx, uint module_degree, uint module_idx){
    FreeModule * module = (FreeModule *) this;

    FreeModuleOperationGeneratorPair operation_generator = module->basis_element_to_opgen_table[module_degree][module_idx];
    uint module_operation_degree = operation_generator.operation_degree;
    uint module_operation_index = operation_generator.operation_degree;
    uint generator_index = operation_generator.generator_index;

    // Now all of the output elements are going to be of the form s * x. Find where such things go in the output vector.
    uint output_block_min = module->generator_to_index_table[module_degree + op_deg][generator_index];
    uint output_block_max = module->generator_to_index_table[module_degree + op_deg][generator_index + 1];
    
    Vector output_block; 
    result->interface->slice(&output_block, result, output_block_min, output_block_max);
    // Now we multiply s * r and write the result to the appropriate position.
    algebra_multiply_basis_elements(module->module.algebra, &output_block, coeff, op_deg, op_idx, module_operation_degree, module_operation_index);
}

void FreeModuleHomomorphism_apply_to_basis_element(FreeModuleHomomorphism * f, Vector * result, uint coeff, uint input_degree, uint input_index){
    VectorInterface * vectImpl = result->interface;
    FreeModuleOperationGeneratorPair operation_generator = f->source->basis_element_to_opgen_table[input_degree][input_index];
    uint operation_degree = operation_generator.operation_degree;
    uint operation_index = operation_generator.operation_degree;
    uint generator_degree = operation_generator.generator_degree;
    uint generator_index = operation_generator.generator_index;

    Vector * output_on_generator = f->outputs[generator_degree][generator_index];
    for(
        VectorIterator it = vectImpl->getIterator(output_on_generator);
        it.has_more; 
        it = vectImpl->stepIterator(it)
    ){
        if(it.value != 0){
            uint c = modPLookup( output_on_generator->p, it.value*coeff);
            module_act_on_basis(f->target, result, c, operation_degree, operation_index, generator_degree, it.index);
        }
    }

}

void FreeModuleAllocateSpaceForNewGenerators(FreeModuleHomomorphism * f, uint num_gens){
    VectorInterface * vectImpl = &f->source->module.algebra->vectorInterface;
    uint p = f->source->module.p;
    uint degree = f->max_computed_degree;
    uint dimension = module_get_dimension(f->target, degree);
    uint vector_size = vectImpl->getSize(p, dimension, 0);
    f->max_computed_degree++;
    f->outputs[degree] = (Vector**)malloc(
        num_gens * sizeof(Vector*) 
        + num_gens * vectImpl->container_size * sizeof(uint64)
        + num_gens * vector_size * sizeof(uint64)
    );
    Vector** vector_ptr_ptr = f->outputs[degree];
    uint64* vector_container_ptr = (uint64*)(vector_ptr_ptr + num_gens);
    uint64* vector_memory_ptr = vector_container_ptr + num_gens* vectImpl->container_size;
    for(uint i = 0; i < num_gens; i++){
        f->outputs[degree][i] = vectImpl->initialize(p, vector_container_ptr, vector_memory_ptr, dimension, 0);
        vector_ptr_ptr ++;
        vector_container_ptr += vectImpl->container_size;
        vector_memory_ptr += vector_size;
    }
}

void addGeneratorToFreeModuleHomomorphism(FreeModuleHomomorphism * f, uint degree, Vector * output){
    assert(f->max_computed_degree == degree);
    assert(output->dimension == module_get_dimension(f->target, degree));
    assert(output->offset == 0);
    VectorInterface * vectImpl = output->interface;
    uint num_gens_in_degree = f->source->number_of_generators_in_degree[degree];
    vectImpl->assign(f->outputs[degree][num_gens_in_degree], output);
    f->source->number_of_generators++;
    f->source->number_of_generators_in_degree[degree]++;
}

// result should be big enough to hold output (how big is that?)
void getHomomorphismMatrix(Vector ** result, FreeModuleHomomorphism * f, uint degree){
    for(int i = 0; i < module_get_dimension(&f->source->module, degree); i++){
        FreeModuleHomomorphism_apply_to_basis_element(f, result[i], 1, degree, i);
    }
}

Kernel * constructKernel(VectorInterface * vectImpl, uint p, uint rows, uint columns){
    Kernel * k = malloc(
        sizeof(Kernel) 
        + columns * sizeof(uint)
        + getMatrixSize(vectImpl, p, rows, columns) * sizeof(uint64)
    );
    k->dimension = rows;
    k->column_to_pivot_row = (uint*)(k + 1);
    k->kernel = initializeMatrix((uint64*)(k->column_to_pivot_row + columns), vectImpl, p, rows, columns);
    return k;
}