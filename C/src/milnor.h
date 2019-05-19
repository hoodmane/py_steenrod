//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H

#include "combinatorics.h"
#include "milnor_datatypes.h"



void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n);
void freeMilnorBasisPpartTable(MilnorAlgebra * algebra);
void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long n );
void freeMilnorBasisQPartTable(MilnorAlgebra * algebra);

void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree);
void freeMilnorBasis(MilnorAlgebra * algebra);
void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);

MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, MonomialIndex idx);
MonomialIndex GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);

void MilnorProductEven(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
void MilnorProductFullQpart(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement m1, unsigned long f);
void MilnorProduct(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r, MilnorBasisElement s);


void initialize_milnor_matrix(unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);
bool step_milnor_matrix(unsigned long  p, unsigned long M[MAX_XI_TAU][MAX_XI_TAU], P_part r, P_part s);

#endif //CSTEENROD_MILNOR_H
