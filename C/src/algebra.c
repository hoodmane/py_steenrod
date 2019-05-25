//
// Created by Hood on 5/20/2019.
//

// TODO: split up this file a bit.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "combinatorics.h"
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
            uint vectorSize = vectorInterface.getSize(p, number_of_generators_in_degree[output_degree]);
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
                    vectorInterface.initialize(p, current_ptr_5, number_of_generators_in_degree[output_degree]);
                    current_ptr_4 ++;
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

void addActionToFiniteDimensionalModule(                                \
    FiniteDimensionalModule * module,                                   \
    uint operation_degree, uint operation_idx,        \
    uint input_degree, uint input_idx,                \
    uint * output                                              \
){
    // (in_deg) -> (out_deg) -> (op_index) -> (in_index) -> Vector
    Vector * output_vector = &module->actions[input_degree][input_degree + operation_degree][operation_idx][input_idx];
    uint p = module->module.p;
    VectorInterface vectorInterface = module->module.algebra->vectorInterface;
    vectorInterface.pack(p, output_vector, output);
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
    uint p = this->p;
    VectorInterface vectorInterface = this->algebra->vectorInterface;
    vectorInterface.add(p, result, &module->actions[mod_degree][mod_degree + op_degree][op_index][mod_index], coeff);
}

/*
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
    uint output_generator_index = module->generator_to_index_table[module_degree + op_deg][generator_index];
    uint output_block_dimension = module->generator_to_index_table[module_degree + op_deg][generator_index + 1] - output_generator_index;
    // Set up a vector pointing into the appropriate block in result that we use multiply_basis_elements to write to.
    Vector output_block;
    output_block.p = result->p;
    output_block.degree = result->degree;
    output_block.dimension = output_block_dimension;
    output_block.vector = result->vector + output_generator_index;
    // Now we multiply s * r and write the result to the appropriate position.
    algebra_multiply_basis_elements(module->module.algebra, &output_block, coeff, op_deg, op_idx, module_operation_degree, module_operation_index);
}

void FreeModuleHomomorphism_apply_to_basis_element(FreeModuleHomomorphism * f, Vector * result, uint coeff, uint input_degree, uint input_index){
    FreeModuleOperationGeneratorPair operation_generator = f->source->basis_element_to_opgen_table[input_degree][input_index];
    uint operation_degree = operation_generator.operation_degree;
    uint operation_index = operation_generator.operation_degree;
    uint generator_degree = operation_generator.generator_degree;
    uint generator_index = operation_generator.generator_index;

    Vector output_on_generator = f->outputs[generator_degree][generator_index];
    for(int i = 0; i < output_on_generator.dimension; i++){
        uint c = output_on_generator.vector[i];
        if(c != 0){
            c = ModPositive(c*coeff, output_on_generator.p);
            module_act_on_basis(f->target, result, c, operation_degree, operation_index, output_on_generator.degree, i);
        }
    }

}

// result should be big enough to hold output (how big is that?)
void getHomomorphismMatrix(uint * Vector, FreeModuleHomomorphism * f, uint degree){
    for(int i = 0; i < module_get_dimension(&f->source->module, degree); i++){
        FreeModuleHomomorphism_apply_to_basis_element(f, result[i], 1, degree, i);
    }
}

#define max(a,b) ((a) > (b) ? (a) : (b))

// Invariants:
//    Kernel needs to be big enough.
void generateOldKernelAndComputeNewKernel(Resolution * resolution, Kernel * kernel, uint degree){
    uint p = resolution->algebra->p;
    uint homological_degree = resolution->internal_degree_to_resolution_stage[degree] + 1;
    FreeModuleHomomorphism current_differential = resolution->resolution_differentials[homological_degree];
    FreeModuleHomomorphism previous_differential = resolution->resolution_differentials[homological_degree - 1];
    FreeModule * source = current_differential.source;
    Module * target = current_differential.target;
    // assert(target == previous_differential.source)
    uint source_dimension = module_get_dimension(&current_differential->source->module, degree);
    uint target_dimension = module_get_dimension(current_differential->target, degree);

    source->max_generator_degree ++;
    // assert(source->max_generator_degree == degree);

    // The Homomorphism matrix has size source_dimension x target_dimension, but we are going to augment it with an
    // identity matrix so that gives a matrix with dimensions source_dimension x (target_dimension + source_dimension).
    // Later we're going to write into this same matrix an isomorphism source/image + new vectors --> kernel
    // This has size target_dimension x (2*target_dimension).
    // This latter matrix may be used to find a preimage of an element under the differential.
    // The information about the kernel is moved to a separate temporary kernel matrix and eventually can be found by looking at
    // the image of the next step in the resolution.
    uint rows = max(source_dimension, target_dimension);
    uint columns = 2*target_dimension + source_dimension;
    Vector * matrix[rows];
    uint matrix_contents[columns*rows];
    memset(matrix_contents, 0, )
    uint * current_ptr;
    for(uint i = 0; i < rows; i++){
        matrix[i]->p = p;
        matrix[i]->degree = -1;
        matrix[i]->dimension = columns;
        matrix[i]->vector = current_ptr;
        current_ptr += columns;
    }
    getHomomorphismMatrix(matrix, current_differential, degree);
    // Write the identity matrix into the right block
    for(int i = 0; i < source_dimension; i++){
        matrix[i][target_dimension + i] = 1;
    }
    // Row reduce
    rows = source_dimension;
    columns = target_dimension + rows;
    uint column_to_pivot_row[columns];
    row_reduce(matrix, column_to_pivot_row, p, rows, columns);

    // Locate first kernel row
    uint first_kernel_row = rows;
    for(uint column = target_dimension; column < target_dimension + source_dimension; column ++){
        if(column_to_pivot_row[column] != -1){
            first_kernel_row = column_to_pivot_row[i];
            break;
        }
    }
    uint kernel_rows = rows - first_kernel_row;
    // Write pivots into kernel
    for(uint column = 0; column < source_dimension; column ++){
        kernel->column_to_pivot_row[column] = column_to_pivot_row[column + target_dimension] - kernel_rows;
    }

    // Copy kernel matrix into kernel
    kernel->dimension = kernel_rows;
    for(uint row = 0; row < kernel_rows; row++){
        memcpy(
            kernel->kernel[row].vector,
            matrix[first_kernel_row + row] + target_dimension,
            source_dimension * sizeof(uint)
        );
    }
    current_differential->kernel[degree] = kernel;

    Vector * previous_kernel = previous_differential->kernel[degree];
    uint current_target_row = first_kernel_row;
    uint kernel_size = 0;
    // Find the quotient previous_kernel/image and add a new generator to source for
    for(uint i = 0; i < target_dimension; i++){
        // Did we find an element of the kernel of the last differential that we haven't hit?
        if(column_to_pivot_row[i] < 0 && previous_differential->kernel_column_to_pivot_row[degree][i] >= 0){
            // Look up the vector that we're missing and add a generator hitting it.
            source->number_of_generators ++;
            source->number_of_generators_in_degree[degree]++;
            Vector kernel_vector = previous_kernel[i];
            memcpy(matrix[current_target_row], kernel_vector.vector, target_dimension * sizeof(uint))
            memset(matrix[current_target_row] + target_dimension, 0, source_dimension * sizeof(uint));
            matrix[current_target_row][target_dimension + source_dimension + kernel_size] = 1;
            current_target_row++;
            kernel_size++;
        }
    }
    // Now the part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate this much space and copy.
    uint coimage_to_image_rows = current_differential;
    uint coimage_to_image_columns = target_dimension + source_dimension + kernel_size;

    // TODO: This allocation code beuints in a separate function. Maybe we can make a single function for allocate and copy.
    Vector* coimage_to_image = malloc(
            coimage_to_image_rows * sizeof(Vector)
            + coimage_to_image_rows * coimage_to_image_columns * sizeof(uint);
    );
    uint * curr_uint_ptr = (uint*)(coimage_to_image + coimage_to_image_rows);
    for(uint row = 0; row < coimage_to_image_rows; row++){
        coimage_to_image[row].p = p;
        coimage_to_image[row].dimension = coimage_to_image_columns;
        coimage_to_image[row].degree = -1; // Does this field mean anything?
        coimage_to_image->vector = curr_uint_ptr;
        curr_uint_ptr += coimage_to_image_columns;
    }
    //assert(curr_uint_ptr == (uint*)(coimage_to_image + coimage_to_image_rows) + coimage_to_image_columns);
    current_differential.coimage_to_image_isomorphism[degree] = coimage_to_image;

    // copy matrix contents to coimage_to_image
    for(uint row = 0; row < coimage_to_image_rows; row++) {
        memcpy(coimage_to_image[i].vector, matrix[i], coimage_to_image_columns * sizeof(uint));
    }


    resolution->internal_degree_to_resolution_stage[degree] ++;
}

//void computeSurjection()


//FreeModuleHomomorphism * constructFreeModuleHomomorphism(FreeModule * source, Module * target);
void addOutputToFreeModuleHomomorphism(FreeModuleHomomorphism * f, uint gen_index, uint * output);


#include "milnor.h"
//
//int main(){
//    MilnorAlgebra * A = constructMilnorAlgebra(2, false, NULL);
//    Algebra * algebra = (Algebra*) A;
//    algebra_compute_basis(algebra, 50);
//    uint max_generator_degree = 4;
//    uint number_of_generators_in_degree[5] = {1,1,1,1,1};
////    constructFiniteDimensionalModule(algebra, max_generator_degree, number_of_generators_in_degree);
//    return 0;
//}

*/