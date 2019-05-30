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

def checkDsquaredZero(source, d_first, d_second):
    for g in source.gens:
        dg = cmodules.c_apply_homomorphism(d_first, source.get_generator(g))
        ddg = cmodules.c_apply_homomorphism(d_second, dg)
        if len(ddg) > 0:
            print(g)

def printMatrixInfo(d_first, d_second, degree):
    M = cmodules.c_homomorphism_to_matrix(d_first, degree)
    for idx, ar in enumerate(M):
        x = cmodules.free_module_elt_from_array(d_first.target, degree, ar)
        print(idx, x)
        if(x==0):
            continue
        dx = cmodules.c_apply_homomorphism(d_second, x)
        if dx != 0:
            print("    ", "dx = ", dx)

if __name__ == "__main__":
    cmilnor.makeCMilnorAlgebra(p=2, degree=50)
    A = steenrod.MilnorAlgebra(p=2)
    Sq = A.Sq
    M = steenrod_module.FiniteSteenrodModule(p=2)
    x0 = M.add_basis_element("x0", 0)
    M.validate()
    cM = cmodules.FiniteSteenrodModule_to_C(M)
    print("resolve")
    res = resolve(M, 20)
    print("done resolving")
    print(res.contents.modules[2].contents.max_generator_degree)
    F0 = cmodules.FreeModule_from_c(res.contents.modules[1], A, "x0")
    F1 = cmodules.FreeModule_from_c(res.contents.modules[2], A, "x1")
    F2 = cmodules.FreeModule_from_c(res.contents.modules[3], A, "x2")
    F3 = cmodules.FreeModule_from_c(res.contents.modules[4], A, "x3")
    F4 = cmodules.FreeModule_from_c(res.contents.modules[5], A, "x4")
    F5 = cmodules.FreeModule_from_c(res.contents.modules[6], A, "x5")
    d1 = cmodules.homomorphism_from_c(res.contents.differentials[2], F1, F0)
    d2 = cmodules.homomorphism_from_c(res.contents.differentials[3], F2, F1)
    d3 = cmodules.homomorphism_from_c(res.contents.differentials[4], F3, F2)
    d4 = cmodules.homomorphism_from_c(res.contents.differentials[5], F4, F3)
    d5 = cmodules.homomorphism_from_c(res.contents.differentials[6], F5, F4)
    x11 = F1.get_generator("x11")
    x12 = F1.get_generator("x12")
    x14 = F1.get_generator("x14")
    x18 = F1.get_generator("x18")
    x116 = F1.get_generator("x116")
    x22 = F2.get_generator("x22")
    x24 = F2.get_generator("x24")
    x25 = F2.get_generator("x25")
    x28 = F2.get_generator("x28")
    x29 = F2.get_generator("x29")
    x210 = F2.get_generator("x210")
    x216 = F2.get_generator("x216")

    for i in range(8):
        op = Sq(i)
        print(op * cmodules.c_apply_homomorphism(d2, x29) - cmodules.c_apply_homomorphism(d2, op * x29))

    #checkDsquaredZero(F2, d2, d1)
    #checkDsquaredZero(F3, d3, d2)
    #printMatrixInfo(d2, d1, 17)

    # The correct cocycle for h_0h_4:
    # Sq(16)*x11  +  Sq(7, 3)*x11  +  Sq(4, 4)*x11  +  Sq(10, 1)*x14  +  Sq(4, 3)*x14  
    # +  Sq(6, 0, 1)*x14  +  Sq(3, 1, 1)*x14  +  Sq(0, 2, 1)*x14  +  Sq(9)*x18 + Sq(1)*x116