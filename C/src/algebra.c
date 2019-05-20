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

void addVector(Vector * target, Vector * source){
    for(long i = 0; i < target->dimension; i++){
        target->vector[i] += source->vector[i];
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
