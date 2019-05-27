//
// Created by Hood on 5/20/2019.
//

// TODO: split up this file a bit.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "combinatorics.h"
#include "algebra.h"
#include "modules.h"
#include "resolution.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

// Invariants:
void generateOldKernelAndComputeNewKernel(Resolution * resolution, uint degree){
    VectorInterface * vectImpl = &resolution->algebra->vectorInterface;
    uint p = resolution->algebra->p;
    uint homological_degree = resolution->internal_degree_to_resolution_stage[degree] + 1;
    FreeModuleHomomorphism * current_differential  = &resolution->resolution_differentials[homological_degree];
    FreeModuleHomomorphism * previous_differential = &resolution->resolution_differentials[homological_degree - 1];
    FreeModule * source = current_differential->source;
    // assert(target == previous_differential.source)
    uint source_dimension = module_getDimension(&current_differential->source->module, degree);
    uint target_dimension = module_getDimension(current_differential->target, degree);

    source->max_generator_degree ++;
    // assert(source->max_generator_degree == degree);

    // The Homomorphism matrix has size source_dimension x target_dimension, but we are going to augment it with an
    // identity matrix so that gives a matrix with dimensions source_dimension x (target_dimension + source_dimension).
    // Later we're going to write into this same matrix an isomorphism source/image + new vectors --> kernel
    // This has size target_dimension x (2*target_dimension).
    // This latter matrix may be used to find a preimage of an element under the differential.
    // The information about the kernel is moved to a separate temporary kernel matrix and eventually can be found by looking at
    // the image of the next step in the resolution.
    uint entries_per_limb = vectImpl->getEntriesPer64Bits(p);
    // Pad the target dimension so that it ends in an aligned position.
    uint padded_target_dimension = 
        ((target_dimension + entries_per_limb - 1)/entries_per_limb)*entries_per_limb;
    uint rows = max(source_dimension, target_dimension);
    uint columns = padded_target_dimension + source_dimension + target_dimension;
    uint64 matrix_memory[getMatrixSize(vectImpl, p, rows, columns)];
    Matrix *full_matrix = initializeMatrix(matrix_memory, vectImpl, p, rows, columns);
    // For the first stage we just want the part of size source_dimension x (padded_target_dimension + source_dimension).
    // Slice matrix out of full_matrix.
    Matrix matrix_contents;
    Matrix *matrix = &matrix_contents;
    matrix->rows = source_dimension;
    matrix->columns = padded_target_dimension + source_dimension; 
    Vector *slice_matrix[rows];
    matrix->matrix = slice_matrix;
    for(uint i = 0; i < matrix->columns; i++){
        vectImpl->slice(slice_matrix[i], full_matrix->matrix[i], 0, matrix->rows);
    }
    FreeModuleHomomorphism_getMatrix(matrix, current_differential, degree);
    // Write the identity matrix into the right block
    for(int i = 0; i < source_dimension; i++){
       vectImpl->setEntry(matrix->matrix[i], target_dimension + i, 1);
    }
    // Row reduce
    rows = source_dimension;
    columns = target_dimension + rows;
    int column_to_pivot_row[columns];
    rowReduce(matrix, column_to_pivot_row);

    // Stage 1: Find kernel of current differential
    // Locate first kernel row
    uint first_kernel_row = rows;
    for(uint i = target_dimension; i < columns; i ++){
        if(column_to_pivot_row[i] != -1){
            first_kernel_row = column_to_pivot_row[i];
            break;
        }
    }
    uint kernel_size = rows - first_kernel_row;
    Kernel * kernel = Kernel_construct(vectImpl, p, kernel_size, target_dimension);
    // Write pivots into kernel
    for(uint i = first_kernel_row; i < rows; i++){
        kernel->column_to_pivot_row[i] = column_to_pivot_row[i + target_dimension] - kernel_size;
    }

    // Copy kernel matrix into kernel
    for(uint row = 0; row < kernel_size; row++){
        Vector slice;
        vectImpl->slice(&slice, matrix->matrix[first_kernel_row + row], padded_target_dimension, padded_target_dimension + source_dimension);
        vectImpl->assign(kernel->kernel->matrix[row], &slice);
    }
    current_differential->kernel[degree] = kernel;

    // Stage 2: Hit kernel of previous differential. 
    Kernel * previous_kernel = previous_differential->kernel[degree];
    // We no longer care about the kernel rows since we stored them somewhere else, 
    // so we're going to write over them.
    uint current_target_row = first_kernel_row;
    // Find basis of quotient previous_kernel/image and add new free module generators to hit
    uint homology_size = 0;
    for(uint i = 0; i < target_dimension; i++){
        if(column_to_pivot_row[i] < 0 && previous_kernel->column_to_pivot_row[i] >= 0){
            // Look up the vector that we're missing and add a generator hitting it.
            Vector * kernel_vector = previous_kernel->kernel->matrix[i];
            Vector slice; 
            vectImpl->slice(&slice, full_matrix->matrix[current_target_row], 0, target_dimension);
            vectImpl->assign(&slice, kernel_vector);

            vectImpl->slice(&slice, full_matrix->matrix[current_target_row], padded_target_dimension, columns);
            vectImpl->setToZero(&slice);        
            vectImpl->setEntry(&slice, source_dimension + homology_size, 1);
            current_target_row++;
            homology_size++;
        }
    }
    free(previous_kernel); // This information can now be found in this differential's coimage_to_image_matrix.
    
    // TODO: make a FreeModule function that increments degree and mallocs space for these vectors.
    FreeModuleHomomorphism_AllocateSpaceForNewGenerators(current_differential, homology_size);
    for(uint i = 0; i < homology_size; i++){
        Vector slice; 
        vectImpl->slice(&slice, full_matrix->matrix[first_kernel_row + i], 0, target_dimension);
        FreeModuleHomomorphism_addGenerator(current_differential, degree,  &slice);
    }
    FreeModule_ConstructBlockOffsetTable(source, degree);

    // Now the part of the matrix that contains interesting information is current_target_row * (target_dimension + source_dimension + kernel_size).
    // Allocate a matrix coimage_to_image with these dimensions.
    uint coimage_to_image_rows = current_target_row;
    uint coimage_to_image_columns = padded_target_dimension + source_dimension + kernel_size;
    Matrix *coimage_to_image = constructMatrix(vectImpl, p, coimage_to_image_rows, coimage_to_image_columns);
    current_differential->coimage_to_image_isomorphism[degree] = coimage_to_image;
    // Copy matrix contents to coimage_to_image
    for(uint i = 0; i < coimage_to_image_rows; i++) {
        Vector slice;
        vectImpl->slice(&slice, full_matrix->matrix[i], 0, coimage_to_image_columns);
        vectImpl->assign(coimage_to_image->matrix[i], &slice);
    }
    int useless_pivot_row_info[coimage_to_image_columns];
    rowReduce(coimage_to_image, useless_pivot_row_info);
    // TODO: assertion about contents of useless_pivot_row_info?
    // Should contain [0,1,2,3,...,n,-1,-1,-1,..., -1].
    resolution->internal_degree_to_resolution_stage[degree] ++;
}


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

