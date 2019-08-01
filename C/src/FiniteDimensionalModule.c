#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "FiniteDimensionalModule.h"

// The allocator is horrendous so we're going to separate it out.
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint name_length, int max_basis_degree, uint *graded_dimension);

FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, char *name, int min_degree, int max_basis_degree, uint *graded_dimension){
    uint name_length = strlen(name) + 1;
    FiniteDimensionalModule *result = FiniteDimensionalModule_allocate(algebra, name_length, max_basis_degree - min_degree, graded_dimension);
    result->module.p = algebra->p;
    strcpy(result->module.name, name);
    result->module.algebra = algebra;
    result->module.computeBasis = FiniteDimensionalModule_computeBasis;
    result->module.getDimension = FiniteDimensionalModule_getDimension;
    result->module.actOnBasis = FiniteDimensionalModule_actOnBasis;
    result->module.min_degree = min_degree;
    result->module.max_degree = INT_MAX; // There is no "max degree" for a finite module -- we've computed it through an infinite range.
    result->module.max_computed_degree = INT_MAX;
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
FiniteDimensionalModule *FiniteDimensionalModule_allocate(Algebra *algebra, uint name_length, int max_basis_degree, uint *graded_dimension){
    uint p = algebra->p;
    uint name_length_padded = ((name_length + sizeof(uint) - 1)/sizeof(uint)) * sizeof(uint);
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
    for(int input_degree = 0; input_degree < max_basis_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            continue;
        }
        for(int output_degree = input_degree + 1; output_degree < max_basis_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                continue;
            }
            uint number_of_operations = Algebra_getDimension(algebra, output_degree - input_degree, input_degree);
            action_matrix_size_3 += sizeof(Vector**) * number_of_operations;
            action_matrix_size_4 += sizeof(Vector*) * number_of_operations * graded_dimension[input_degree];
            uint vector_size = Vector_getSize(p, graded_dimension[output_degree], 0);
            action_matrix_size_5 += vector_size * number_of_operations * graded_dimension[input_degree];
        }
    }

    action_matrix_size_2 += action_matrix_size_1;
    action_matrix_size_3 += action_matrix_size_2;
    action_matrix_size_4 += action_matrix_size_3;
    action_matrix_size_5 += action_matrix_size_4;

    FiniteDimensionalModule *result = malloc(
            sizeof(FiniteDimensionalModule)
            + name_length_padded
            + max_basis_degree * sizeof(uint) // graded_dimension
            + action_matrix_size_5
    );
    result->module.name = (char *) (result + 1);
    result->graded_dimension = (uint *) (result->module.name + name_length_padded);
    char *top_of_action_table = (char *)(result->graded_dimension + max_basis_degree);
    Vector *****current_ptr_1 = (Vector *****) top_of_action_table;
    Vector ****current_ptr_2 = (Vector ****) (top_of_action_table + action_matrix_size_1);
    Vector ***current_ptr_3 = (Vector ***) (top_of_action_table + action_matrix_size_2);
    Vector **current_ptr_4 = (Vector **) (top_of_action_table + action_matrix_size_3);
    char *current_ptr_5 = (char *) (top_of_action_table + action_matrix_size_4);
    result->actions = current_ptr_1; 
    memset(top_of_action_table, 0, action_matrix_size_4);
    for(int input_degree = 0; input_degree < max_basis_degree; input_degree++){
        if(graded_dimension[input_degree] == 0){
            current_ptr_1 ++;
            continue;
        }
        *current_ptr_1 = current_ptr_2;
        current_ptr_2 += input_degree + 1;
        for(int output_degree = input_degree + 1; output_degree < max_basis_degree; output_degree++){
            if(graded_dimension[output_degree] == 0){
                current_ptr_2 ++;
                continue;
            }
            *current_ptr_2 = current_ptr_3;
            uint number_of_operations = Algebra_getDimension(algebra, output_degree - input_degree, input_degree);
            for(uint operation_idx = 0; operation_idx < number_of_operations; operation_idx ++){
                *current_ptr_3 = current_ptr_4;
                for(uint input_idx = 0; input_idx < graded_dimension[input_degree]; input_idx ++ ){
                    *current_ptr_4 = Vector_initialize(p, &current_ptr_5, graded_dimension[output_degree], 0);
                    current_ptr_4 ++;
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
    int operation_degree, uint operation_index,
    int input_degree, uint input_index,
    uint *output
){
    input_degree -= module->module.min_degree;
    assert(input_degree >= 0);
    uint output_degree = input_degree + operation_degree;
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    // printf("operation_degree: %d, module_degree: %d, output_degree: %d, operation_index: %d, module_index: %d\n", operation_degree, input_degree, output_degree, operation_index, input_index);    
    Vector *output_vector = module->actions[input_degree][output_degree][operation_index][input_index];
    Vector_pack(output_vector, output);
}


Vector *FiniteDimensionalModule_getAction(
    FiniteDimensionalModule *module,
    int operation_degree, uint operation_index,
    int module_degree, uint module_index
){
    module_degree -= module->module.min_degree;
    int output_degree = module_degree + operation_degree;
    // if(operation_degree == 4 && operation_index == 0){
        // printf(">>operation_degree: %d, module_degree: %d, output_degree: %d, operation_index: %d, module_index: %d\n", operation_degree, module_degree, output_degree, operation_index, module_index);    
    // }    
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
    FiniteDimensionalModule *module = (FiniteDimensionalModule*)this;
    assert(op_index < Algebra_getDimension(this->algebra, op_degree, mod_degree));
    assert(mod_index < Module_getDimension(this, mod_degree));
    uint output_dimension = Module_getDimension(this, mod_degree + op_degree);    
    if(mod_degree + op_degree >= module->max_basis_degree || output_dimension == 0){
        return;
    }  
    // Why do we take this slice?
    char output_block_memory[Vector_getSize(this->p, 0, 0)];    
    char *output_block_ptr = output_block_memory;
    Vector *output_block = Vector_initialize(this->p, &output_block_ptr, 0, 0);     
    Vector_slice(output_block, result, 0, output_dimension); 
    Vector *output = FiniteDimensionalModule_getAction(module, op_degree, op_index, mod_degree, mod_index);
    // if(op_degree == 4 && op_index == 0 && mod_degree - this->min_degree == 0 && mod_index == 0){
        // printf("    fdm_aob -- op_degree: %d, op_index: %d, mod_degree: %d, mod_index: %d\n", op_degree, op_index, mod_degree, mod_index);        
        // Vector_print("    ---- %s\n", output);
    // }
    Vector_add(output_block, output, coeff);
}