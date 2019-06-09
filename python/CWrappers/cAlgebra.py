from cMilnorAlgebra import *
from cAdemAlgebra import *

c_algebra_constructors = {}
c_algebra_constructors["MilnorAlgebra"] = cMilnorAlgebra
c_algebra_constructors["AdemAlgebra"] = cAdemAlgebra

algebra_instances = {}

def getAlgebra(algebra_type, *, p, generic=None, max_degree=0):
    if generic == None:
        generic = p!=2
    key = (algebra_type, p, generic)
    if key in algebra_instances:
        algebra = algebra_instances[key]
    else:
        algebra = c_algebra_constructors[algebra_type](p, generic)
        algebra_instances[key] = algebra
    algebra.generateBasis(max_degree)
    return algebra


if __name__ == "__main__":
    A = getAlgebra("AdemAlgebra", p=3)
    A.generateBasis(20)
    P = A.P
    A.multiply(P(1),P(1))        