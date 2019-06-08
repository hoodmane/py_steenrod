import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes import *
from ctypes_wrap import *
from cFpVector import *

import steenrod

class cAdemBasisElement:
    def __init__(self, algebra, elt):
        self.c_elt = elt
        self.algebra = algebra
        self.generateReprString(algebra.generic)
        self.degree = self.c_elt.contents.degree
        assert self.degree < algebra.max_degree
        self.idx = CSteenrod.AdemAlgebra_basisElement_toIndex(self.algebra.c_algebra, self.c_elt.contents)
        self.py_rep = None

    def generateReprString(self, generic):
        P_or_Sq = "P" if generic else "Sq"
        ps = []
        bocksteins = self.c_elt.contents.bocksteins
        for i in range(self.c_elt.contents.P_length):
            if self.algebra.generic and (bocksteins & 1) !=0:
                ps += ["b"]
            bocksteins >>= 1
            ps += [P_or_Sq + str(self.c_elt.contents.Ps[i])]
        if self.algebra.generic and (bocksteins & 1) != 0:
            ps += ["b"]
        repr_str = " ".join(ps)
        if(len(repr_str) == 0):
            repr_str = "0"
        self.repr_str = repr_str
        
    def __repr__(self):
        return self.repr_str

    def mult(self, other):
        if(self.algebra != other.algebra):
            raise Exception("Algebras don't match.")
        if(self.degree + other.degree > self.algebra.max_degree):
            raise Exception("Degree exceeds maximum.")
        result = cAdemElement(algebra = self.algebra, degree = self.degree + other.degree)
        CSteenrod.AdemAlgebra_multiply(
            self.algebra.c_alg_ptr, 
            result.cVect.vector, 1, self.degree, self.idx, other.degree, other.idx, 0
        )
        return result

    def toPyTuple(self):
        if self.py_rep:
            return self.py_rep
        length = self.c_elt.contents.P_length
        if self.algebra.generic:
            bocksteins = [0] * (length + 1)
            bitstring = self.c_elt.contents.bocksteins
            while bitstring != 0:
                if bitstring & 1 != 0:
                    bocksteins[i] = 1
                bitstring >>= 1
        Ps = [0] * length
        for i in range(length):
            Ps[i] = self.c_elt.contents.Ps[i]
        if self.algebra.generic:
            result = [0] * (2* length + 1)
            result[0::2] = bocksteins
            result[1::2] = Ps
        else:
            result = tuple(Ps)
        self.py_rep = result
        return result

    def toPy(self):
        return self.algebra.py_algebra.get_basis_element(self.toPyTuple())

class cAdemElement:
    def __init__(self, *, algebra, degree, vector = None):
        self.algebra = algebra
        self.degree = degree
        self.dimension = algebra.getDimension(degree)
        if(vector):
            if vector.dimension != self.dimension:
                raise Exception("Dimensions don't match!")
            self.cVect = vector
        else:
            self.cVect = cVector(self.algebra.p, dim = self.dimension)
        self.repr_string = None
        self.freed = False

    def generateReprString(self):
        basis_strs = []
        basis = self.algebra.getBasis(self.degree)
        for idx, v in enumerate(self.cVect.unpack()):
            if v==0:
                continue
            cur_str = ""                
            if v!=1:
                cur_str += "%s*" % v
            cur_str += str(basis[idx])
            basis_strs += [cur_str]
        self.repr_string = " + ".join(basis_strs)
        if len(self.repr_string) == 0:
            self.repr_string = "0" 


    def __repr__(self):
        if not self.repr_string:
            self.generateReprString()
        return self.repr_string

    def add(self, elt, c=1):
        self.repr_string = None
        if elt.degree != self.degree:
            raise Error("Degrees don't match.")
        c = c % self.algebra.p
        self.cVect.add(elt.cVect, c)
    
    def __setitem__(self, idx, value):
        self.repr_string = None
        self.cVect[idx] = value

    def __getitem__(self, idx, value):
        return self.cVect[idx] 

    def __iter__(self):
        return self.cVect.__iter__()

    def toPy(self):
        result = {}
        basis = self.algebra.getBasis(self.degree)
        for i, entry in enumerate(self.cVect.unpack()):
            if entry == 0:
                continue
            result[basis[i].toPyTuple()] = entry
        return self.algebra.py_algebra.get_element(result)        

    def free(self):
        self.freed = True
        self.cVect.free()

class cAdemAlgebra:
    def __init__(self, p, generic=None, max_degree=0):
        if generic is None:
            generic = p != 2
        self.p = p
        self.generic = generic
        self.c_algebra =  CSteenrod.AdemAlgebra_construct(p, generic, POINTER(c_Profile)())
        self.py_algebra = steenrod.AdemAlgebra.getInstance(p, generic)
        self.c_alg_ptr = cast(self.c_algebra, POINTER(c_Algebra))
        self.max_degree = 0
        self.generateBasis(max_degree)
        self.bases = {}
        self.Sq = self.py_algebra.Sq
        self.P = self.py_algebra.P
        self.b = self.py_algebra.b

    def generateBasis(self, max_degree):
        if(max_degree > self.max_degree):
            CSteenrod.AdemAlgebra_generateBasis(self.c_alg_ptr, max_degree)
            self.max_degree = max_degree

    def getDimension(self, degree):
        return CSteenrod.AdemAlgebra_getDimension(self.c_alg_ptr, degree, 0)

    def getBasis(self, degree):
        if degree in self.bases:
            return self.bases[degree]
        if degree > self.max_degree:
            raise Exception("C basis only known through degree %s < %s." % (self.max_degree, degree))
        c_basisElementList = CSteenrod.AdemAlgebra_getBasis(self.c_alg_ptr, degree, 0)
        result = [None] * int(c_basisElementList.length)
        for i in range(c_basisElementList.length):
            b = c_basisElementList.list[i]
            result[i] = cAdemBasisElement(self, b)
        self.bases[degree] = result
        return result
    
    def basis_elt_to_elt(self, b):
        idx = CSteenrod.AdemAlgebra_basisElement_toIndex(self.c_algebra, b.c_elt)
        v = cVector(self.p, dim=self.getDimension(b.degree))
        v[idx] = 1
        return cAdemElement(algebra=self,degree=b.degree, vector=v)

    def multiply(self, m1, m2):
        m1_deg = m1.degree()
        m2_deg = m2.degree()
        out_degree = m1_deg + m2_deg
        if out_degree > self.max_degree:
            raise Exception("C basis only known through degree %s < %s." % (self.max_degree, out_degree))
        result = cAdemElement(algebra = self, degree=out_degree)
        b1 = next(iter(m1))
        b2 = next(iter(m2))
        m1_idx = self.idxFromPy(b1)
        m2_idx = self.idxFromPy(b2)
        CSteenrod.AdemAlgebra_multiply(self.c_alg_ptr, result.cVect.vector, 1, m1_deg, m1_idx, m2_deg, m2_idx,0)
        py_res = result.toPy()
        # result.free()
        return py_res

    def idxFromPy(self, b):
        bitstring = 0
        p_part = b
        if self.generic:
            for idx, x in enumerate(b[0::2]):
                bitstring |= x << idx
            p_part = b[1::2]
        p_array = (c_uint * len(p_part))(*p_part)
        basis_elt = c_AdemBasisElement()
        basis_elt.degree = self.py_algebra.basis_degree(b)
        basis_elt.bocksteins = bitstring
        basis_elt.P_length = len(p_part)
        basis_elt.Ps = p_array
        basis_elt_ptr = POINTER(c_AdemBasisElement)(basis_elt)
        return CSteenrod.AdemAlgebra_basisElement_toIndex(self.c_algebra, basis_elt_ptr)
        
    def basisEltFromPy(self, b):
        degree = self.py_algebra.basis_q_degree(b) + self.py_algebra.basis_p_degree(b)
        idx = self.idxFromPy(b)
        return self.getBasis(degree)[idx]
                







if __name__ == "__main__":
    pass
    A = cAdemAlgebra(p=2, max_degree=100)
    Sq = A.py_algebra.Sq
    P = A.py_algebra.P
    #A3 = makeCAdemAlgebra(p=3, dim=200)
    
