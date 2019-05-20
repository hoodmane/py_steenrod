//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H
#include <stdbool.h>

#include "combinatorics.h"
#include "milnor_datatypes.h"
#include "algebra.h"



void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n);
void freeMilnorBasisPpartTable(MilnorAlgebra * algebra);
void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long n );
void freeMilnorBasisQPartTable(MilnorAlgebra * algebra);

void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree);
void freeMilnorBasis(MilnorAlgebra * algebra);
void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);

MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, unsigned long degree, unsigned long index);
unsigned long GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);

void MilnorProductEven(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
void MilnorProductFullQpart(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement m1, unsigned long f);
void MilnorProduct(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement r, MilnorBasisElement s);


void initialize_milnor_matrix(unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);
bool step_milnor_matrix(unsigned long  p, unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);

Algebra * getMilnorAlgebra(unsigned long p, bool generic, Profile profile);
#endif //CSTEENROD_MILNOR_H
