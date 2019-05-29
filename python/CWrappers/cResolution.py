import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cmodules
import cFpVector
import cmilnor
import steenrod
import steenrod_module

def resolve(module, degree):
    cRes = CSteenrod.Resolution_construct(module.c_module, degree, degree)
    CSteenrod.Resolution_resolveThroughDegree(cRes, degree)
    return cRes

if __name__ == "__main__":
    cmilnor.makeCMilnorAlgebra(p=2, degree=50)
    A = steenrod.MilnorAlgebra(p=2)
    Sq = A.Sq
    M = steenrod_module.FiniteSteenrodModule(p=2)
    x0 = M.add_basis_element("x0", 0)
    M.validate()
    cM = cmodules.FiniteSteenrodModule_to_C(M)
    print("resolve")
    res = resolve(M, 10)
    print("done resolving")
    F0 = cmodules.FreeModule_from_c(res.contents.modules[1], A, "x0")
    F1 = cmodules.FreeModule_from_c(res.contents.modules[2], A, "x1")
    F2 = cmodules.FreeModule_from_c(res.contents.modules[3], A, "x2")
    F3 = cmodules.FreeModule_from_c(res.contents.modules[4], A, "x3")
    d1 = cmodules.homomorphism_from_c(res.contents.differentials[2], F1, F0)
    d2 = cmodules.homomorphism_from_c(res.contents.differentials[3], F2, F1)
    d3 = cmodules.homomorphism_from_c(res.contents.differentials[4], F3, F2)
    x11 = F1.get_generator("x11")
    x12 = F1.get_generator("x12")
    x14 = F1.get_generator("x14")
    x22 = F2.get_generator("x22")
    x24 = F2.get_generator("x24")
    x25 = F2.get_generator("x25")