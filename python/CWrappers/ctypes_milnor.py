from ctypes import *
from ctypes_milnor_datatypes import *

def wrap_milnor(CSteenrod):
    #void generateMilnorBasisPpartTable(MilnorAlgebra * algebra, unsigned long n);
    #void generateMilnorBasisQpartTable(MilnorAlgebra * algebra, unsigned long n );
    
    #void GenerateMilnorBasis2(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
    #void GenerateMilnorBasisGeneric(MilnorAlgebra * algebra, unsigned long old_n, unsigned long n);
    #void GenerateMilnorBasis(MilnorAlgebra * algebra, unsigned long max_degree);
    CSteenrod.GenerateMilnorBasis.argtypes = [POINTER(c_MilnorAlgebra), c_ulong]
    
    #MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, unsigned long degree, unsigned long idx);
    CSteenrod.GetMilnorBasisElementFromIndex.argtypes = [POINTER(c_MilnorAlgebra), c_ulong, c_ulong]
    CSteenrod.GetMilnorBasisElementFromIndex.restype = c_MilnorBasisElement
    #unsigned long GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);
    CSteenrod.GetIndexFromMilnorBasisElement.argtypes = [POINTER(c_MilnorAlgebra), c_MilnorBasisElement]
    CSteenrod.GetIndexFromMilnorBasisElement.restype = c_ulong   
    
    #void MilnorProductEven(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement r_elt, MilnorBasisElement s_elt);
    #void MilnorProductFullQpart(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement m1, unsigned long f);
    #void MilnorProduct(MilnorAlgebra * algebra, Vector * result, MilnorBasisElement r, MilnorBasisElement s);
    CSteenrod.MilnorProduct.argtypes = [POINTER(c_MilnorAlgebra), POINTER(c_Vector), c_MilnorBasisElement, c_MilnorBasisElement]
    
    #unsigned long** initialize_milnor_matrix(P_part r, P_part s);
    #bool step_milnor_matrix(unsigned long  p, unsigned long ** M, P_part r, P_part s);
