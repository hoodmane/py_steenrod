//
// Created by Hood on 5/8/2019.
//

#ifndef CSTEENROD_MILNOR_H
#define CSTEENROD_MILNOR_H

#include "milnor_datatypes.h"



void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n);
void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long n );


void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree);
void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);

MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, MonomialIndex idx);
MonomialIndex GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);

MilnorElement * MilnorProductEven(MilnorAlgebra * algebra, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
MilnorElement * MilnorProductFullQpart(MilnorAlgebra * algebra , MilnorBasisElement m1, unsigned long f);
MilnorElement *  MilnorProduct(MilnorAlgebra * algebra, MilnorBasisElement r, MilnorBasisElement s);

unsigned long** initialize_milnor_matrix(P_part r, P_part s);
bool step_milnor_matrix(unsigned long  p, unsigned long ** M, P_part r, P_part s);

#endif //CSTEENROD_MILNOR_H
