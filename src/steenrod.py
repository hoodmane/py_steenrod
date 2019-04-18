import itertools
from functools import reduce

from memoized import memoized
import adem
import milnor
from FpVectorSpace import Vector

def implementedByAssignmentLaterInThisFile():
    assert False, "We implement this by assignment from Vector.linearly_extend_map later in this file Steenrod.py"


class MilnorElement(Vector):
    def __init__(self, algebra, dictionary):
        self.algebra = algebra
        super(MilnorElement, self).__init__(algebra.p, dictionary)
        
    def __mul__(self, v):
        implementedByAssignmentLaterInThisFile()
            
    def to_adem(self):
        implementedByAssignmentLaterInThisFile()
            
    def basis_elt_to_string(self, b):
        if(self.algebra.generic):
            Qs = ["Q(%s)" % s for s in b[0]]
            Ps = ""
            if len(b[1]) > 0:
                Ps = "P(%s)" % ", ".join([str(s) for s in b[1]])
            return " ".join(Qs + [Ps]) or '1'
        else:
            Sqs = ""
            if len(b) > 0:
                Sqs = "Sq(%s)" % ", ".join([str(s) for s in b])
            return Sqs or '1'
            
MilnorElement.__mul__ = Vector.linearly_extend_map(milnor.product)
       
       
class AdemElement(Vector):
    def __init__(self, dictionary, *,  algebra):
        self.algebra = algebra
        super(AdemElement, self).__init__(algebra.p, dictionary)
    
    def __mul__(self, v):
        implementedByAssignmentLaterInThisFile()
       
    def to_milnor(self):
        implementedByAssignmentLaterInThisFile()
               
    def basis_elt_to_string(self, b):
        if(self.algebra.generic):
            Ps = b[1::2]
            bs = b[0::2]
            Ps = [ " P%s" % s for s in Ps]
            bs = [ " b" if s else "" for s in b[0::2]]
            result = [None]*len(b)
            result[1::2] = Ps
            result[0::2] = bs
            return "".join(result)[1:]
        else:
            return " ".join(["Sq%s" % s for s in b])

AdemElement.__mul__ = Vector.linearly_extend_map(adem.product)


class AdemAlgebra:
    def __init__(self, p, generic = None):
        self.p = p
        self.generic = generic
    
    @staticmethod
    def getAlgebra(p, generic = None):
        if generic is None:
            generic = p != 2
        if (p,generic) not in AdemAlgebra.instance_dict:
            AdemAlgebra(self, p, generic)
        return AdemAlgebra.instance_dict[(p,generic)]
    
    def zero(self):
        return AdemElement(self, {})

    def unit():
        return AdemElement(self, {() : 1})    
    
    def b():
        return AdemElement(self, {(1,) : 1})
    
    def Sq(n):
        if self.generic:
            return AdemElement(self, {(n % 2, n//2, 0) : 1})
        else:
            return AdemElement(self, {(n,) : 1})
    
    def P(n):
        if self.generic:
            return AdemElement(self, {(0, n, 0) : 1})
        else:
            return AdemElement(self, {(2*n,) : 1})
    
    def bP(n):
        if self.generic:
            return AdemElement(self, {(1, n, 0) : 1})
        else:
            return AdemElement(self, {(2*n + 1,) : 1})
        
AdemAlgebra.instance_dict = {}
    
    
class MilnorAlgebra:
    def __init__(self, p, generic = None, profile = None):
        if generic is None:
            generic = p != 2
        self.p = p
        self.generic = generic
        self.profile = profile
        self.unit_monomial = ((),()) if generic else ()
    
    @staticmethod
    def getAlgebra(p, generic = None, profile = None, truncation = None):
        if generic is None:
            generic = p != 2
        instance_dict = MilnorAlgebra.instance_dict
        if (p,generic, profile) not in instance_dict:
            instance_dict[(p,generic, profile)] = MilnorAlgebra(p, generic, profile)
        return instance_dict[(p,generic, profile)]
    
    def zero(self):
        return MilnorElement(self, {})
    
    def unit(self):
        return MilnorElement(self, { self.unit_monomial : 1 })
    
    def Sq(self, *l):
        if self.generic:
            raise NotImplementedError()
        else:
            return MilnorElement(self, { l : 1 })
    
    def Q(self, *l):
        if self.generic:
            return MilnorElement(self, { (l,()) : 1 })
        else:
            return reduce(AdemElement.__mul__, [MilnorElement(self, { ((0,)* (i) + (1,)) : 1 }) for i in l])
    
    def P(self, *l):
        if self.generic:
            return MilnorElement(self, {((), l) : 1})
        else:
            return MilnorElement(self, {tuple(2 * i for i in l) : 1})

MilnorAlgebra.instance_dict = {}



#@memoized
def adem_to_milnor_on_basis_2(b):
    A = MilnorAlgebra.getAlgebra(2)
    unit = A.unit()
    sqs = [A.Sq(i) for i in b]
    return reduce(lambda a,b : a*b, sqs, unit)

#@memoized
def adem_to_milnor_on_basis_generic(b, p):
    A = MilnorAlgebra.getAlgebra(p)
    unit = A.unit()
    Ps = b[1::2]
    bs = b[0::2]
    Q0 = A.Q(0)
    Ps = [A.P(j) for j in Ps]
    bs = [Q0 if j else None for j in bs]
    sqs = [None]*len(b)
    sqs[1::2] = Ps
    sqs[0::2] = bs
    sqs = [j for j in result if j is not None]
    return reduce(lambda a,b : a*b, sqs, unit)
    
def adem_to_milnor_on_basis(b, *, p, generic):
    if generic:
        return adem_to_milnor_on_basis_generic(b, p)
    else:
        return adem_to_milnor_on_basis_2(b)

#@memoized
def milnor_to_adem_on_basis_2(b):
    t = [0] * len(b)
    t[-1] = b[-1]
    for i in range(len(b) - 2, -1, -1):
        t[i] = b[i] + 2 * t[i + 1]
    t = tuple(t)
    x = adem_to_milnor_on_basis_2(t)
    x.pop(b)
    result = x.to_adem()
    result[t] = 1
    return result

#@memoized
def milnor_to_adem_on_basis_generic(b, *, p):
    e = b[0]
    s = b[1]
    pad_length = max(*e)
    s += (0,) * (pad_length - len(s)) 
    t = [0,] * (2*len + 1)
    
    for i in e:
        t[2*i] = 1
            
    t[-2] = s[-1] + t[-1]
    for i in range(2, s.length + 1):
        t[-2*i] = t[-2*i + 1] + s[-i] + p * t[-2*i + 2]
    t = tuple(t)
    x = adem_to_milnor_on_basis_generic(t, p)
    x.pop(b)
    x.scale_in_place(-1)
    result = x.to_adem()
    result[t] = 1
    return result

def milnor_to_adem_on_basis(b, *, p, generic):
    if generic:
        return milnor_to_adem_on_basis_generic(b, p)
    else:
        return milnor_to_adem_on_basis_2(b)

AdemElement.to_milnor = Vector.linearly_extend_map(adem_to_milnor_on_basis)
MilnorElement.to_adem = Vector.linearly_extend_map(milnor_to_adem_on_basis)

