import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes_wrap import *
import cAlgebra
import cFiniteDimensionalModule
import cFreeModule
import cFreeModuleHomomorphism
import cFpVector
import steenrod
import steenrod_module

def resolve(module, degree):
    cRes = CSteenrod.Resolution_construct(module.c_module, degree, 0, 0)
    CSteenrod.Resolution_resolveThroughDegree(cRes, degree)
    return cRes

def printDimensionsToFile(cRes, filename):
    buffer = (c_char * 1000000)()
    CSteenrod.Resolution_gradedDimensionString(buffer, cRes)
    charp = c_char_p(buffer[:])
    s = charp.value.decode("utf-8")
    steenrod_module.write_file(filename, s)


def checkDsquaredZero(source, d_first, d_second):
    for g in source.gens:
        dg = cFreeModuleHomomorphism.apply(d_first, source.get_generator(g))
        ddg = cFreeModuleHomomorphism.apply(d_second, dg)
        if len(ddg) > 0:
            print(g)

def printMatrixInfo(d_first, d_second, degree):
    M = cFreeModuleHomomorphism.toMatrix(d_first, degree)
    for idx, ar in enumerate(M):
        x = FreeModule.elementfromArray(d_first.target, degree, ar)
        print(idx, x)
        if(x==0):
            continue
        dx = cFreeModuleHomomorphism.apply(d_second, x)
        if dx != 0:
            print("    ", "dx = ", dx)

def checkDifferential(d, x):
    for (op, g) in x:
        prod1 = op* d(d.source.get_generator(g))
        prod2 = d(d.source.get_basis_element(op, g))
        if len(prod1 - prod2) != 0:
            print("op:", op, "g:", g, ":")
            print("    ", prod1)
            print("    ", prod2)
            print("    ", prod1  + prod2)
            return prod1  + prod2

if __name__ == "__main__":
    degree = 75
    p = 2
    algebra_type = "AdemAlgebra"
    filename = "%s%s_%s" % (algebra_type, p, degree)
    A = cAlgebra.getAlgebra(algebra_type, p=p, max_degree=degree)
    # P = A.P
    # bP = A.py_algebra.bP
    # b = A.py_algebra.b()
    M = steenrod_module.FiniteSteenrodModule(p=p)
    x0 = M.add_basis_element("x0", 0)
    M.validate()
    cM = cFiniteDimensionalModule.toC(M, algebra_type)
    res = resolve(M, degree)
    # print("done resolving")
    # # print(res.contents.modules[2].contents.max_generator_degree)
    # res_modules = []
    # for i in range(degree-1):
    #     F = cFreeModule.fromC(res.contents.modules[i+1], A, "x" + str(i) + "_")
    #     res_modules.append(F)
    #     globals()["F" + str(i)] = F
    #     for i in F.gens:
    #         globals()[i] = F.get_generator(i)

    # for i in range(1, degree-1):
    #     globals()["d" + str(i)] = cFreeModuleHomomorphism.fromC(res.contents.differentials[i+1], res_modules[i], res_modules[i-1])

    #printDimensionsToFile(res, "output_data/%s.json" % filename)


    # for i in range(8):
    #     op = Sq(i)
    #     print(op * cmodules.c_apply_homomorphism(d2, x29) - cmodules.c_apply_homomorphism(d2, op * x29))

    #checkDsquaredZero(F2, d2, d1)
    #checkDsquaredZero(F3, d3, d2)
    #printMatrixInfo(d2, d1, 17)

    # The correct cocycle for h_0h_4:
    # Sq(16)*x11  +  Sq(7, 3)*x11  +  Sq(4, 4)*x11  +  Sq(10, 1)*x14  +  Sq(4, 3)*x14  
    # +  Sq(6, 0, 1)*x14  +  Sq(3, 1, 1)*x14  +  Sq(0, 2, 1)*x14  +  Sq(9)*x18 + Sq(1)*x116