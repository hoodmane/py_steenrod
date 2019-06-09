import os,sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from ctypes import *
from ctypes_wrap import *
from cFpVector import *
import steenrod

class cMilnorBasisElement:
    def __init__(self, algebra, elt):
        self.c_elt = elt
        self.algebra = algebra
        self.generateReprString(algebra.generic)
        self.degree = self.c_elt.p_degree + self.c_elt.q_degree
        self.idx = CSteenrod.MilnorAlgebra_basisElement_toIndex(self.algebra.c_algebra, self.c_elt)
        self.py_rep = None

    def generateReprString(self, generic):
        repr_list = []
        if generic and self.c_elt.q_part != 0:
            repr_list += ["Q(%s)" % self.c_elt.q_part]
        if self.c_elt.p_length != 0:
            ps = [None] * self.c_elt.p_length
            for i in range(self.c_elt.p_length):
                ps[i] = str(self.c_elt.p_part[i])
            pstr = ", ".join(ps)
            pstr = "%s(%s)" % ("P" if generic else "Sq",  pstr)
            repr_list += [pstr]
        repr_str = "*".join(repr_list)
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
        result = cMilnorElement(algebra = self.algebra, degree = self.degree + other.degree)
        CSteenrod.MilnorAlgebra_multiply(
            self.algebra.c_alg_ptr, 
            result.cVect.vector, 1, self.degree, self.idx, other.degree, other.idx, 0
        )
        return result

    def toPyTuple(self):
        if self.py_rep:
            return self.py_rep
        q_part = ()
        i = 0
        bitstring = self.c_elt.q_part
        while bitstring != 0:
            if bitstring & 1 != 0:
                q_part += (i,)
            i += 1
            bitstring >>= 1
        p_part = [0] * self.c_elt.p_length
        for i in range(self.c_elt.p_length):
            p_part[i] = self.c_elt.p_part[i]
        if self.algebra.generic:
            result = (q_part, tuple(p_part))
        else:
            result = tuple(p_part)
        self.py_rep = result
        return result

    def toPy(self):
        return self.algebra.py_algebra.get_basis_element(self.toPyTuple())

class cMilnorElement:
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

class cMilnorAlgebra:
    def __init__(self, p, generic=None, max_degree=None):
        if generic is None:
            generic = p != 2
        self.p = p
        self.generic = generic
        if (p, generic) in cMilnorAlgebra.cAlgebras:
            self.c_algebra = cMilnorAlgebra.cAlgebras[(p, generic)]
        else:
            self.c_algebra =  CSteenrod.MilnorAlgebra_construct(p, generic, POINTER(c_Profile)())
            cMilnorAlgebra.cAlgebras[(p, generic)] = self.c_algebra
        self.py_algebra = steenrod.MilnorAlgebra.getInstance(p, generic)
        self.c_alg_ptr = cast(self.c_algebra, POINTER(c_Algebra))
        self.basis_type = cMilnorBasisElement
        self.max_degree = 0
        self.bases = {}
        self.Sq = self.py_algebra.Sq
        self.P = self.py_algebra.P
        self.Q = self.py_algebra.Q
        if max_degree:
            self.generateBasis(max_degree)        

    def generateBasis(self, max_degree):
        if(max_degree > self.max_degree):
            print("py generating basis:", max_degree)
            CSteenrod.MilnorAlgebra_generateBasis(self.c_alg_ptr, max_degree)
            self.max_degree = max_degree

    def getDimension(self, degree):
        return CSteenrod.MilnorAlgebra_getDimension(self.c_alg_ptr, degree, 0)

    def getBasis(self, degree):
        if degree in self.bases:
            return self.bases[degree]
        if degree > self.max_degree:
            raise Exception("C basis only known through degree %s < %s." % (self.max_degree, degree))
        c_basisElementList = CSteenrod.MilnorAlgebra_getBasis(self.c_alg_ptr, degree, 0)
        result = [None] * int(c_basisElementList.length)
        for i in range(c_basisElementList.length):
            b = c_basisElementList.list[i]
            result[i] = cMilnorBasisElement(self, b)
        self.bases[degree] = result
        return result
    
    def basis_elt_to_elt(self, b):
        idx = CSteenrod.MilnorAlgebra_basisElement_toIndex(self.c_algebra, b.c_elt)
        v = cVector(self.p, dim=self.getDimension(b.degree))
        v[idx] = 1
        return cMilnorElement(algebra=self,degree=b.degree, vector=v)

    def multiply(self, m1, m2):
        m1_deg = m1.degree()
        m2_deg = m2.degree()
        out_degree = m1_deg + m2_deg
        if out_degree > self.max_degree:
            raise Exception("C basis only known through degree %s < %s." % (self.max_degree, out_degree))
        result = cMilnorElement(algebra = self, degree=out_degree)
        b1 = next(iter(m1))
        b2 = next(iter(m2))
        m1_idx = self.idxFromPy(b1)
        m2_idx = self.idxFromPy(b2)
        CSteenrod.MilnorAlgebra_multiply(self.c_alg_ptr, result.cVect.vector, 1, m1_deg, m1_idx, m2_deg, m2_idx,0)
        py_res = result.toPy()
        result.free()
        return py_res

    def basisEltFromIdx(self, degree, idx):
        return self.getBasis(degree)[idx]        

    def idxFromPy(self, b):
        bitstring = 0
        p_part = b
        if self.generic:
            for x in b[0]:
                bitstring |= 1 << x
            p_part = b[1]
        p_array = (c_uint * len(p_part))(*p_part)
        mbe = c_MilnorBasisElement()
        mbe.q_degree = self.py_algebra.basis_q_degree(b)
        mbe.q_part = bitstring
        mbe.p_degree = self.py_algebra.basis_p_degree(b)
        mbe.p_length = len(p_part)
        mbe.p_part = p_array
        degree = mbe.q_degree + mbe.p_degree
        return CSteenrod.MilnorAlgebra_basisElement_toIndex(self.c_algebra, mbe)
        
    def basisEltFromPy(self, b):
        degree = self.py_algebra.basis_q_degree(b) + self.py_algebra.basis_p_degree(b)
        idx = self.idxFromPy(b)
        return self.getBasis(degree)[idx]

cMilnorAlgebra.cAlgebras = {}







if __name__ == "__main__":
    pass
    A = cMilnorAlgebra(p=2, max_degree=100)
    Sq = A.py_algebra.Sq
    #A3 = makeCMilnorAlgebra(p=3, dim=200)
    
