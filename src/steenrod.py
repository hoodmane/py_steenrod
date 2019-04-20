import itertools
from functools import reduce

from memoized import memoized
import adem
import milnor
from FpVectorSpace import Vector

def implementedByAssignmentLaterInThisFile():
    assert False, "We implement this by assignment from Vector.linearly_extend_map later in this file Steenrod.py"


class AdemElement(Vector):
    def __init__(self, dictionary, *,  algebra):
        self.algebra = algebra
        self.module = algebra
        super(AdemElement, self).__init__(algebra.p, dictionary)
    
    def __mul__(self, v):
        if type(v) == AdemElement:
            return self.multiply(v)
        elif callable(getattr(v, "adem_act", None)):
            return v.adem_act(self)
        elif type(v) == int:
            result = self.algebra.getElement(self)
            result.scale_in_place(v)
            return result
        else:
            raise TypeError()
    
    def __rmul__(self, v):
#        if type(v) == AdemElement:
#            return self.multiply(v)
#        elif callable(getattr(v, "adem_act", None)):
#            return v.adem_act(self)
#        elif type(v) == int:
        if type(v) == int:
            result = self.algebra.getElement(self)
            result.scale_in_place(v)
            return result
        else:
            raise TypeError()
    
    def multiply(self, v):
        implementedByAssignmentLaterInThisFile()
       
    def to_milnor(self):
        implementedByAssignmentLaterInThisFile()
        
    def basis_degree(self, b):
        if self.generic:
            result  = 2*(self.p - 1) * sum(b[1::2])
            result += sum(b[::2])
        else:
            result = sum(b)
        return result
        
    def basis_elt_to_string(self, basis_elt):
        if(self.algebra.generic):
            result = adem.adem_basis_elt_generic_map(P_fn = lambda P : " P%s" % P, b = " b", basis_elt = basis_elt)
            return "".join(result)[1:]
        else:
            return " ".join(["Sq%s" % s for s in basis_elt]) # This is adem.adem_2_map inlined

AdemElement.multiply = Vector.linearly_extend_map(adem.product, kw_param = "algebra")

class MilnorElement(Vector):
    def __init__(self, algebra, dictionary):
        self.algebra = algebra
        self.module = algebra
        super(MilnorElement, self).__init__(algebra.p, dictionary)
        
    def __mul__(self, v):
        if type(v) == MilnorElement:
            return self.multiply(v)
        elif callable(getattr(v, "milnor_act", None)):
            return v.milnor_act(self)
        elif type(v) == int:
            result = self.algebra.getElement(self)
            result.scale_in_place(v)
            return result
        else:
            raise TypeError()
            
    def __rmul__(self, v):
#        if type(v) == AdemElement:
#            return self.multiply(v)
#        elif callable(getattr(v, "adem_act", None)):
#            return v.adem_act(self)
#        elif type(v) == int:
        if type(v) == int:
            result = self.algebra.getElement(self)
            result.scale_in_place(v)
            return result
        else:
            raise TypeError()
                        
    def basis_degree(self, b):
        raise NotImplementedError()
    
    def multiply(self, v):        
        implementedByAssignmentLaterInThisFile()
            
    def to_adem(self):
        implementedByAssignmentLaterInThisFile()
            
    def basis_elt_to_string(self, basis_elt):
        if(self.algebra.generic):
            Qs = ["Q(%s)" % s for s in basis_elt[0]]
            Ps = ""
            if len(basis_elt[1]) > 0:
                Ps = "P(%s)" % ", ".join([str(s) for s in basis_elt[1]])
            return " ".join(Qs + [Ps]) or '1'
        else:
            Sqs = ""
            if len(basis_elt) > 0:
                Sqs = "Sq(%s)" % ", ".join([str(s) for s in basis_elt])
            return Sqs or '1'
            
MilnorElement.multiply = Vector.linearly_extend_map(milnor.product, kw_param = "algebra")
       

class AdemAlgebra:
    def __init__(self, p, generic = None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic
    
    @staticmethod
    def getInstance(p, generic = None):
        if generic is None:
            generic = p != 2
        if (p,generic) not in AdemAlgebra.instance_dict:
            AdemAlgebra.instance_dict[(p,generic)] = AdemAlgebra(p, generic)
        return AdemAlgebra.instance_dict[(p,generic)]
    
    @staticmethod
    def getInstanceFromAlgebra(algebra):
        return AdemAlgebra.getInstance(algebra.p, algebra.generic)
    
    def getElement(self, d):
        return AdemElement(d, module = self)
        
    def getBasisElement(self, b):
        return AdemElement({b : 1}, algebra = self)
    
    def basis(self, n):
        adem.basis(n, algebra = self)
    
    def zero(self):
        return AdemElement({}, algebra = self)

    def unit(self):
        return AdemElement({() : 1}, algebra = self)  
    
    def b(self):
        return AdemElement({(1,) : 1}, algebra = self)
        
    def b_or_unit(self, epsilon):
        if self.generic and epsilon:
            return self.b()
        else:
            return self.unit()
    
    def Sq(self, n):
        if self.generic:
            return AdemElement({(n % 2, n//2, 0) : 1}, algebra = self)
        else:
            return AdemElement({(n,) : 1}, algebra = self)
    
    def P(self, n):
        if self.generic:
            return AdemElement({(0, n, 0) : 1}, algebra = self)
        else:
            return AdemElement({(2*n,) : 1}, algebra = self)
    
    def bP(self, n):
        if self.generic:
            return AdemElement({(1, n, 0) : 1}, algebra = self)
        else:
            return AdemElement({(2*n + 1,) : 1}, algebra = self)
        
AdemAlgebra.instance_dict = {}
    
    
class MilnorAlgebra(milnor.MinimalMilnorAlgebra):
    def __init__(self, p, generic = None, profile = None, truncation = None):
        super(MilnorAlgebra, self).__init__(p, generic, profile, truncation)
        print(str(self.generic))
        self.unit_monomial = ((),()) if generic else ()
    
    @staticmethod
    def getInstance(p, generic = None, profile = None, truncation = None):
        if generic is None:
            generic = p != 2
        instance_dict = MilnorAlgebra.instance_dict
        if (p,generic, profile) not in instance_dict:
            instance_dict[(p,generic, profile)] = MilnorAlgebra(p, generic, profile, truncation)
        return instance_dict[(p,generic, profile)]
    
    def getInstanceFromAlgebra(algebra):
        return MilnorAlgebra.getInstance(algebra.p, algebra.generic, algebra.getattr("profile", None), algebra.getattr("truncation", None))
    
    def getElement(self, d):
        return MilnorElement(d, algebra = self)

    def getBasisElement(self, b):
        return MilnorElement({b : 1}, algebra = self)        
        
    def basis(self, n):
        return [MilnorElement({b:1}, algebra = self) for b in milnor.basis(n, algebra = self)]
    
    def zero(self):
        return MilnorElement({}, algebra = self)
    
    def unit(self):
        return MilnorElement({ self.unit_monomial : 1 }, algebra = self)
    
    def Sq(self, *l):
        if self.generic:
            raise NotImplementedError()
        else:
            return MilnorElement({ l : 1 }, algebra = self)
    
    def Q(self, *l):
        if self.generic:
            return MilnorElement({ (l,()) : 1 }, algebra = self)
        else:
            return reduce(AdemElement.__mul__, [MilnorElement( { ((0,)* (i) + (1,)) : 1 }, algebra = self) for i in l])
    
    def P(self, *l):
        if self.generic:
            return MilnorElement({((), l) : 1}, algebra = self)
        else:
            return MilnorElement({tuple(2 * i for i in l) : 1}, algebra = self)

MilnorAlgebra.instance_dict = {}


#@memoized
def adem_to_milnor_on_basis_2(b):
    A = MilnorAlgebra.getInstance(2)
    unit = A.unit()
    sqs = [A.Sq(i) for i in b]
    return reduce(lambda a,b : a*b, sqs, unit)

#@memoized
def adem_to_milnor_on_basis_generic(basis_elt, p):
    A = MilnorAlgebra.getInstance(p)
    unit = A.unit()
    sqs = adem.adem.adem_basis_elt_generic_map(P_fn = A.P, b = Q0, basis_elt = basis_elt)
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

def milnor_to_adem_on_basis(b, *, algebra):
    if algebra.generic:
        return milnor_to_adem_on_basis_generic(b, p)
    else:
        return milnor_to_adem_on_basis_2(b)

AdemElement.to_milnor = Vector.linearly_extend_map(adem_to_milnor_on_basis, kw_param = "algebra")
MilnorElement.to_adem = Vector.linearly_extend_map(milnor_to_adem_on_basis, kw_param = "algebra")



def adem_antipode_on_basis(basis_elt, *, algebra):
    antipode = algebra.unit()
    milnor_alg = MilnorAlgebra.getInstance(algebra.p, algebra.generic)
    if not algebra.generic:
        for n in basis_elt:
            antipode = self(sum(milnor_alg.basis(n))) * antipode
    else:
        adem.adem_basis_elt_generic_map(P_fn = lambda P : P, b = -algebra.Q(0), basis_elt = basis_elt)
        # ... not sure what happens here
        Ps = b[1::2]
        bs = b[0::2]
        [-algebra.Q(0) for i in bs if i != 0]
    
        B = milnor_alg.basis(n * 2 * (p-1))
        s = self(0)
        for b in B:
            if len(b.leading_support()[0]) == 0:
                s += self(b)
        antipode = (-1)**n * s * antipode
    return antipode
