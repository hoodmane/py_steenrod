import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cFpVector
import cAlgebra
import cFreeModule
import steenrod
import steenrod_module
from FreeModule import *


def toC(f, c_S, c_T):
    c_f = CSteenrod.FreeModuleHomomorphism_construct(c_S, cast(c_T, POINTER(c_Module)), 20)
    for i in range(c_S.contents.max_generator_degree+1):
        CSteenrod.FreeModuleHomomorphism_AllocateSpaceForNewGenerators(
            f, i, c_S.contents.number_of_basis_elements_in_degree[i]
        )
    
    for (g, deg) in c_S.py_module.gens.items():
        CSteenrod.FreeModule_getDimension(sour)
        c_output = cFpVector.vector_to_C()
        CSteenrod.FreeModule_setOutput(f, g, c_S.generator_indices[g], )

class cHomomorphism:
    def __init__(self, *, cf, source, target):
        self.cf = cf
        self.source = source
        self.target = target
    
    def __call__(self, input):
        return self.apply(input)

    def toMatrix(f, degree):
        c_source_cast = cast(f.source.c_module, POINTER(c_Module))
        c_target_cast = cast(f.target.c_module, POINTER(c_Module))
        input_dimension = CSteenrod.FreeModule_getDimension(c_source_cast, degree)
        output_dimension = CSteenrod.FreeModule_getDimension(c_target_cast, degree)
        M = cFpVector.cMatrix(2, input_dimension, output_dimension)
        CSteenrod.FreeModuleHomomorphism_getMatrix(f.cf, M.cM, degree)
        result = M.unpack()
        M.free()
        return result
        
    def apply(f, element):
        degree = element.degree()
        c_source = f.source.c_module
        c_target_cast = cast(f.target.c_module, POINTER(c_Module))
        output_dimension = CSteenrod.FreeModule_getDimension(c_target_cast, degree)
        c_result = cFpVector.cVector(f.source.p, output_dimension)
        for ((elt_op, gen), coeff) in element.items():
            elt_op_deg = elt_op.degree()
            elt_op = next(iter(elt_op))
            elt_op_idx = f.source.c_algebra.idxFromPy(elt_op)
            gen_deg = f.source.gens[gen]
            gen_idx = f.source.generator_indices[gen]
            # print(elt_op_deg, elt_op_idx, gen_deg, gen_idx)
            elt_idx = CSteenrod.FreeModule_operationGeneratorToIndex(f.source.c_module, elt_op_deg, elt_op_idx, gen_deg, gen_idx)
            # print(elt_op, gen)
            # print("degree: ", degree, "idx: ", elt_idx)
            CSteenrod.FreeModuleHomomorphism_applyToBasisElement(f.cf, c_result.vector, coeff, degree, elt_idx)
            # print(free_module_elt_from_c(f.target, degree, c_result))
        result = cFreeModule.elementFromC(f.target, degree, c_result.vector)
        c_result.free()
        return result

def fromC(cf, S, T):
    f = cHomomorphism(cf=cf, source=S, target=T)
    return f

if __name__ == "__main__":
    A = steenrod.MilnorAlgebra(p=2)
    Sq = A.Sq
    M = steenrod_module.FiniteSteenrodModule(p=2)
    x0 = M.add_basis_element("x0", 0)
    x1 = M.add_basis_element("x1", 1)
    x2 = M.add_basis_element("x2", 2)
    x3 = M.add_basis_element("x3", 3)
    x4 = M.add_basis_element("x4", 4)
    M.add_Sq_action(2, x0, x2)
    M.add_Sq_action(2, "x2", x4)
    M.add_Sq_action(1, "x0", x1)
    M.add_Sq_action(2, "x1", x3)
    M.add_Sq_action(1, "x3", x4)
    M.add_Sq_action(3, "x1", x4)
    M.validate()
    cM = FiniteSteenrodModule_to_C(M)
    c_act_on_fdmodule(cM, Sq(3), x0) # 0
    c_act_on_fdmodule(cM, Sq(0,1), x0) # x3

    T = FreeModule(algebra=A)
    x00 = T.add_generator("x00", 0)
    S = FreeModule(algebra=A)
    x11 = S.add_generator("x11", 1)
    x12 = S.add_generator("x12", 2)
    x14 = S.add_generator("x14", 4)
    x18 = S.add_generator("x18", 8)
    cS = FreeModule_to_c(S, 20)

    d1 = ModuleHomomorphism(S, T)
    d1.add_value("x11", A.Sq(1)*x00)
    d1.add_value("x12", A.Sq(2)*x00)
    d1.add_value("x14", A.Sq(4)*x00)
    d1.add_value("x18", A.Sq(8)*x00)
    # A = F.milnor_algebra