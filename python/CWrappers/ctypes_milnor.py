from ctypes import *
from ctypes_milnor_datatypes import *

def wrap_milnor(CSteenrod):
    #void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n);
    #void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long n );
    
    #void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
    #void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
    #void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree);
    CSteenrod.GenerateMilnorBasis.argtypes = [POINTER(c_MilnorAlgebra), c_ulong]
    
    # MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, MonomialIndex idx);
    CSteenrod.GetMilnorBasisElementFromIndex.argtypes = [POINTER(c_MilnorAlgebra), c_MonomialIndex]
    CSteenrod.GetMilnorBasisElementFromIndex.restype = c_MilnorBasisElement
    # MonomialIndex GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);
    CSteenrod.GetIndexFromMilnorBasisElement.argtypes = [POINTER(c_MilnorAlgebra), c_MilnorBasisElement]
    CSteenrod.GetIndexFromMilnorBasisElement.restype = c_MonomialIndex    
    
    #void MilnorProductEven(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
    #void MilnorProductFullQpart(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement m1, unsigned long f);
    #void MilnorProduct(MilnorAlgebra * algebra, MilnorElement * result, MilnorBasisElement r, MilnorBasisElement s);
    CSteenrod.MilnorProduct.argtypes = [POINTER(c_MilnorAlgebra), POINTER(c_MilnorElement), c_MilnorBasisElement, c_MilnorBasisElement]
    
    #unsigned long** initialize_milnor_matrix(P_part r, P_part s);
    #bool step_milnor_matrix(unsigned long  p, unsigned long ** M, P_part r, P_part s);
