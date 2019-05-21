//
// Created by Hood on 5/20/2019.
//

#include "combinatorics.h"
#include "algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


Vector * allocateVector(unsigned long p, unsigned long degree, unsigned long dimension){
    Vector * result = (Vector*)malloc(sizeof(Vector) + dimension * sizeof(long));
    result->p = p;
    result->degree = degree;
    result->dimension = dimension;
    result->vector = (long*)(result + 1);
    memset(result->vector, 0, dimension * sizeof(long));
    return result;
}

void freeVector(Vector * vector){
    free(vector);
}

void addBasisElementToVector(Vector * target, unsigned long idx, long coeff){
    target->vector[idx] += coeff;
    target->vector[idx] = ModPositive(target->vector[idx], target->p);
}

void addVectors(Vector * target, Vector * source){
    for(long i = 0; i < target->dimension; i++){
        target->vector[i] = ModPositive(source->vector[i] + target->vector[i], target->p);
    }
}

void addArrayToVector(Vector * target, long * array){
    for(long i = 0; i < target->dimension; i++){
        target->vector[i] = ModPositive(array[i] + target->vector[i], target->p);
    }
}

void addVectorToArray(long * target, Vector * source){
    for(long i = 0; i < source->dimension; i++){
        target[i] = ModPositive(target[i] + source->vector[i], source->p);
    }
}

void assignVector(Vector * target, Vector * source){
    memcpy(target->vector, source->vector, target->dimension * sizeof(long));
}
void scaleVector(Vector * target, long s){
    for(long i = 0; i < target->dimension; i++){
        target->vector[i] *= s;
        target->vector[i] = ModPositive(target->vector[i], target->p);
    }
}


bool FiniteDimensionalModule_compute_basis(Module *this, unsigned long dimension);
unsigned long FiniteDimensionalModule_get_basis_dimension(Module* this, unsigned long degree);
void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index);

// The allocator is horrendous so we're going to separate it out.
FiniteDimensionalModule * allocateFiniteDimensionalModule(Algebra * algebra, unsigned long max_generator_degree, unsigned long * number_of_generators_in_degree);

FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, unsigned long max_generator_degree, unsigned long * number_of_generators_in_degree){
    max_generator_degree ++;
    FiniteDimensionalModule * result = allocateFiniteDimensionalModule(algebra, max_generator_degree, number_of_generators_in_degree);
    result->module.p = algebra->p;
    result->module.algebra = algebra;
    result->module.compute_basis = FiniteDimensionalModule_compute_basis;
    result->module.get_basis_dimension = FiniteDimensionalModule_get_basis_dimension;
    result->module.act_on_basis = FiniteDimensionalModule_act_on_basis;
    result->max_degree = max_generator_degree;
    result->dimension = 0;
    for(long i = 0; i < max_generator_degree; i++){
        result->dimension += number_of_generators_in_degree[i];
    }
    memcpy(result->number_of_basis_elements_in_degree, number_of_generators_in_degree, max_generator_degree * sizeof(unsigned long));\
    return result;
}

void freeFiniteDimensionalModule(FiniteDimensionalModule * module){
    free(module);
}

// This is the grossest allocator.
FiniteDimensionalModule * allocateFiniteDimensionalModule(Algebra * algebra, unsigned long max_generator_degree, unsigned long * number_of_generators_in_degree){
    // Count number of triples (x, y, op) with |x| + |op| = |y|.
    // The amount of memory we need to allocate is:
    // # of input_degrees  * sizeof(***Vector)
    // + # of nonempty input degrees * # of output degrees * sizeof(**Vector)
    // + Sum over (nonempty in_deg < nonempty out_deg) of (
    //              # of operations in (out_deg - in_deg) * sizeof(*Vector)
    //              # of operations in (out_deg - in_deg) * # of gens in degree in_degree * sizeof(Vector)
    //              # of operations in (out_deg - in_deg) * # of gens in degree in_degree * # of gens in degree out_degree * sizeof(unsigned long)
    // )
    size_t action_matrix_size_1 = 0, action_matrix_size_2 = 0,
            action_matrix_size_3 = 0, action_matrix_size_4 = 0,
            action_matrix_size_5 = 0;
    unsigned long number_of_nonempty_degrees = 0;
    for(long degree = 0; degree < max_generator_degree; degree++) {
        if(number_of_generators_in_degree[degree] != 0){
            number_of_nonempty_degrees ++;
        }
    }
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> (out_index) -> value
    //  ****    -> ***       -> **Vector   -> *Vector    -> Vector -> long
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
            unsigned long number_of_operations = algebra_get_basis_dimension(algebra, output_degree - input_degree);
            action_matrix_size_3 += sizeof(Vector*)       * number_of_operations;
            action_matrix_size_4 += sizeof(Vector)        * number_of_operations * number_of_generators_in_degree[input_degree];
            action_matrix_size_5 += sizeof(unsigned long) * number_of_operations * number_of_generators_in_degree[input_degree] * number_of_generators_in_degree[output_degree];
        }
    }

    action_matrix_size_2 += action_matrix_size_1;
    action_matrix_size_3 += action_matrix_size_2;
    action_matrix_size_4 += action_matrix_size_3;
    action_matrix_size_5 += action_matrix_size_4;

    FiniteDimensionalModule * result = malloc(
            sizeof(FiniteDimensionalModule)
            + max_generator_degree * sizeof(unsigned long)
            + action_matrix_size_5
    );
    result->number_of_basis_elements_in_degree = (unsigned long * ) (result + 1);
    void * top_of_action_table = (void *) (result + 1) + max_generator_degree * sizeof(unsigned long);
    Vector **** current_ptr_1 = (Vector ****) top_of_action_table;
    Vector *** current_ptr_2 = (Vector ***) (top_of_action_table + action_matrix_size_1);
    Vector ** current_ptr_3 = (Vector **) (top_of_action_table + action_matrix_size_2);
    Vector * current_ptr_4 = (Vector *) (top_of_action_table + action_matrix_size_3);
    long * current_ptr_5 = (long *) (top_of_action_table + action_matrix_size_4);
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
            unsigned long number_of_operations = algebra_get_basis_dimension(algebra, output_degree - input_degree);
            for(int operation_idx = 0; operation_idx < number_of_operations; operation_idx ++){
                *current_ptr_3 = current_ptr_4;
                for(int input_idx = 0; input_idx < number_of_generators_in_degree[input_degree]; input_idx ++ ){
                    current_ptr_4->vector = current_ptr_5;
                    current_ptr_4->p = algebra->p;
                    current_ptr_4->degree = output_degree;
                    current_ptr_4->dimension = number_of_generators_in_degree[output_degree];
                    current_ptr_4 ++;
                    current_ptr_5 += number_of_generators_in_degree[output_degree];
                }
                current_ptr_3 ++;
            }
            current_ptr_2 ++;
        }
        current_ptr_1 ++;
    }
//    assert((long)current_ptr_1 == (long)(top_of_action_table + action_matrix_size_1));
//    assert((long)current_ptr_2 == (long)(top_of_action_table + action_matrix_size_2));
//    assert((long)current_ptr_3 == (long)(top_of_action_table + action_matrix_size_3));
//    assert((long)current_ptr_4 == (long)(top_of_action_table + action_matrix_size_4));
//    assert((long)current_ptr_5 == (long)(top_of_action_table + action_matrix_size_5));
    return result;
}

void addActionToFiniteDimensionalModule(                                \
    FiniteDimensionalModule * module,                                   \
    unsigned long operation_degree, unsigned long operation_idx,        \
    unsigned long input_degree, unsigned long input_idx,                \
    unsigned long * output                                              \
){
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector output_vector = module->actions[input_degree][input_degree + operation_degree][operation_idx][input_idx];
    memcpy(output_vector.vector, output, output_vector.dimension * sizeof(long));
}

bool FiniteDimensionalModule_compute_basis(Module *this, unsigned long dimension){
    return true;
}

unsigned long FiniteDimensionalModule_get_basis_dimension(Module* this, unsigned long degree){
    FiniteDimensionalModule * module = (FiniteDimensionalModule *) this;
    if(degree <= module->max_degree ){
        return module->number_of_basis_elements_in_degree[degree];
    }
    return 0;
}

void FiniteDimensionalModule_act_on_basis(Module * this, Vector *result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index){
    FiniteDimensionalModule* module = ((FiniteDimensionalModule*)this);
    addVectors(result, &module->actions[mod_degree][mod_degree + op_degree][op_index][mod_index]);
}

bool FreeModule_compute_basis(Module* this, unsigned long degree);
unsigned long FreeModule_get_basis_dimension(Module* this, unsigned long degree);
void FreeModule_act_on_basis(Module * this, Vector * result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_idx);

FreeModule * initializeFreeModule(FreeModule * module, Algebra * algebra){
    module->module.algebra = algebra;
    module->module.p = algebra->p;
    module->module.compute_basis = FreeModule_compute_basis;
    module->module.get_basis_dimension = FreeModule_get_basis_dimension;
    module->module.act_on_basis = FreeModule_act_on_basis;
    return module;
}

bool FreeModule_compute_basis(Module* this, unsigned long degree){
    return true;
}

unsigned long FreeModule_get_basis_dimension(Module* this, unsigned long degree){
    FreeModule * module = (FreeModule*) this;
    unsigned long result = 0;
    for(int i = 0; i <= degree; i++){
        result += module->number_of_generators_in_degree[i] * algebra_get_basis_dimension(this->algebra, degree - i);
    }
    return result;
}


void FreeModule_act_on_basis(Module * this, Vector * result, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_idx){

}


//FreeModuleHomomorphism * constructFreeModuleHomomorphism(FreeModule * source, Module * target);
void addOutputToFreeModuleHomomorphism(FreeModuleHomomorphism * f, unsigned long gen_index, long * output);


#include "milnor.h"

int main(){
    MilnorAlgebra * A = constructMilnorAlgebra(2, false, NULL);
    Algebra * algebra = (Algebra*) A;
    algebra_compute_basis(algebra, 50);
    unsigned long max_generator_degree = 4;
    unsigned long number_of_generators_in_degree[5] = {1,1,1,1,1};
//    constructFiniteDimensionalModule(algebra, max_generator_degree, number_of_generators_in_degree);
    return 0;
}