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
        self.adem_element = True
        super(AdemElement, self).__init__(algebra.p, dictionary)
    
    def adem_act(self,v):
        if type(v) != AdemElement:
            raise TypeError()
        if v.algebra.p != self.algebra.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()
        return v.multiply(self)
    
    def multiply(self, v):
        implementedByAssignmentLaterInThisFile()
       
    def to_milnor(self):
        implementedByAssignmentLaterInThisFile()
        
    def basis_degree(self, b):
        if self.module.generic:
            result  = 2*(self.p - 1) * sum(b[1::2])
            result += sum(b[::2])
        else:
            result = sum(b)
        return result
        
    def basis_elt_to_string(self, basis_elt):
        if(self.algebra.generic):
            result = adem.adem_basis_elt_generic_map(P_fn = lambda P : " P%s" % P, b = " b", basis_elt = basis_elt)
            result = "".join(result)[1:]
        else:
            result =  " ".join(["Sq%s" % s for s in basis_elt]) # This is adem.adem_2_map inlined
        return result or "1"
AdemElement.multiply = Vector.linearly_extend_map(adem.product, kw_param = "algebra")

class MilnorElement(Vector):
    def __init__(self, dictionary, *, algebra):
        self.algebra = algebra
        self.module = algebra
        self.milnor_element = True
        super(MilnorElement, self).__init__(algebra.p, dictionary)
        
    def milnor_act(self,v):
        if type(v) != MilnorElement:
            raise TypeError()
        if v.p != self.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()    
        return v.multiply(self)
                        
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
        return AdemElement(d, algebra = self)
        
    def getBasisElement(self, b):
        return AdemElement({b : 1}, algebra = self)
    
    def getInadmissiblePairs(self, max_degree):
        P = self.P if generic else self.Sq
        for relation_dim in range(2, max_degree):
            # We want Sqi*Sqj inadmissible so that means i < 2 * j.
            # relation_dim = i + j so j = relation_dim - i
            # so i < 2 * (relation_dim - i) so i < 2 * relation_dim / 3. 
            # We need to round up so that Python includes the last integer 
            # if 2 * relation_dim / 3 is not an integer
            for i in range(1, int(math.ceil((p * relation_dim) / (p+1)))):
                j = relation_dim - i
                yield (relation_dim, P([i]), P([relation_dim - i]))
            if generic:
                for i in range(1, int(math.ceil((p * relation_dim + 1) / (p+1)))):
                    j = relation_dim - i - 1
                    yield (relation_dim, P([i]), self.bP([relation_dim - i]))                
    
    def basis(self, n):
        return [AdemElement({b:1}, algebra = self) for b in adem.basis(n, algebra = self)]
    
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
    
    def Sq(self, *n):
        if self.generic:
            raise TypeError("Only use A.Sq if p = 2 and A is not generic.")
        # Sq(0) is the unit
        if len(n) == 1 and n[0] == 0:
            n = ()
        return AdemElement({n : 1}, algebra = self)
    
    def P(self, *n):
        if not self.generic:
            raise TypeError("Only use A.P if A is generic.")
        # P(0) is the unit
        if len(n) == 1 and n[0] == 0:
            n = []        
        l = [0] * (2*len(n) + 1)
        l[1::2] = n            
        return AdemElement({tuple(l) : 1}, algebra = self)        
    
    def bP(self, *n):
        if not self.generic:
            raise TypeError("Only use A.bP if A is generic.")
        l = [0] * (2*len(n) + 1)
        l[1::2] = n
        l[0] = 1
        return AdemElement({tuple(l) : 1}, algebra = self)   
    
AdemAlgebra.instance_dict = {}
    
    
class MilnorAlgebra(milnor.MinimalMilnorAlgebra):
    def __init__(self, p, generic = None, profile = None, truncation = None):
        super(MilnorAlgebra, self).__init__(p, generic, profile, truncation)
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

    def basis_even(self, n):
        return [MilnorElement({((), b):1}, algebra = self) for b in milnor.basis_even(n, self.p, self.profile)]
    
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
    sqs = adem.adem_basis_elt_generic_map(P_fn = A.P, b = A.Q(0), basis_elt = basis_elt)
    return reduce(lambda a,b : a*b, sqs, unit)
    
def adem_to_milnor_on_basis(b, *, algebra):
    if algebra.generic:
        return adem_to_milnor_on_basis_generic(b, algebra.p)
    else:
        return adem_to_milnor_on_basis_2(b)

#@memoized
def milnor_to_adem_on_basis_2(b):
    adem_alg = AdemAlgebra.getInstance( p = 2, generic = False )
    if len(b) == 0:
        return adem_alg.unit()
    t = [0] * len(b)
    t[-1] = b[-1]
    for i in range(len(b) - 2, -1, -1):
        t[i] = b[i] + 2 * t[i + 1]
    t = tuple(t)
    x = adem_to_milnor_on_basis_2(t)
    x.pop(b)
    result = x.to_adem(output_module = adem_alg)
    result[t] = 1
    return result

#@memoized
def milnor_to_adem_on_basis_generic(b, *, p):
    (e, s) = b
    adem_alg = AdemAlgebra.getInstance( p = p, generic = True )
    if len(e) == 0 and len(s) == 0:
        return adem_alg.unit()
    if e  == (0,) and len(s) == 0:
        return adem_alg.b()
    pad_length = max(e) if e else 0
    s += (0,) * (pad_length - len(s)) 
    t = [0,] * (2*len(s) + 1)
    
    for i in e:
        t[2*i] = 1
            
    t[-2] = s[-1] + t[-1]
    for i in range(2, len(s) + 1):
        t[-2*i] = t[-2*i + 1] + s[-i] + p * t[-2*i + 2]
    t = tuple(t)
    x = adem_to_milnor_on_basis_generic(t, p)
    x.pop(b)
    x.scale_in_place(-1)
    result = x.to_adem(output_module = adem_alg)
    result[t] = 1
    return result

def milnor_to_adem_on_basis(b, *, algebra):
    if algebra.generic:
        return milnor_to_adem_on_basis_generic(b, p = algebra.p)
    else:
        return milnor_to_adem_on_basis_2(b)

AdemElement.to_milnor = Vector.linearly_extend_map(adem_to_milnor_on_basis, kw_param = "algebra", 
                                get_output_module = lambda alg : MilnorAlgebra.getInstance(alg.p, generic = alg.generic))
MilnorElement.to_adem = Vector.linearly_extend_map(milnor_to_adem_on_basis, kw_param = "algebra",
                                get_output_module = lambda alg : AdemAlgebra.getInstance(alg.p, generic = alg.generic))

def adem_antipode_on_basis(basis_elt, *, algebra):
    antipode = algebra.unit()
    p = algebra.p
    generic = algebra.generic
    milnor_alg = MilnorAlgebra.getInstance(p, generic)
    if not generic:
        elts = [ Vector.sum(milnor_alg.basis(n)) for n in basis_elt ]
    else:
        def P_fn(P):
            return (-1)**P * Vector.sum( milnor_alg.basis_even( P ) )
        
        elts = adem.adem_basis_elt_generic_map(
            P_fn = P_fn,
            b = (-1) * milnor_alg.Q(0), 
            basis_elt = basis_elt
        )
    # Multiply in the reverse order because antipode is an antihomomorphism
    milnor_antipode = reduce( lambda a, b: b * a, elts, milnor_alg.unit() )
    return milnor_antipode.to_adem()

AdemElement.antipode = Vector.linearly_extend_map( adem_antipode_on_basis, kw_param = "algebra")


if __name__ == "__main__":
    A   = AdemAlgebra(p = 2)
    A3  = AdemAlgebra(p = 3)
    Am  = MilnorAlgebra(p = 2)
    A3m = MilnorAlgebra(p = 3) 
    
    print(type(A.Sq(2).to_milnor) == MilnorElement) 
    print(type(Am.Sq(2).to_adem) == AdemElement)
    print(type(A3.P(2).to_milnor) == MilnorElement) 
    print(type(A3m.P(2).to_adem) == AdemElement)
