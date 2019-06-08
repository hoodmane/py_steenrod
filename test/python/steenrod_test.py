import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),os.pardir,"python"))

from steenrod import *


Adem2 = AdemAlgebra(p = 2)
Milnor2 = MilnorAlgebra(p = 2)

Adem3 = AdemAlgebra(p = 3)
Milnor3 = MilnorAlgebra(p = 3)

Sq = Adem2.Sq
P = Adem3.P
bP = Adem3.bP
beta = Adem3.b()

def test_equals():
    assert Sq(3) + Sq(2)*Sq(1) == Sq(2)*Sq(1) + Sq(3)

def test_adem_mult():
    assert type(Sq(0) * Sq(5)) == AdemElement
    assert Sq(0) * Sq(5) == Sq(5)
    assert Sq(1) * Sq(1) == Adem2.zero()
    assert Sq(1) * Sq(2) == Sq(3)
    assert Sq(2) * Sq(1) == Sq(2,1)
    assert Sq(2) * Sq(2) == Sq(3,1)
    assert Sq(3) * Sq(2) == Adem2.zero()
    assert Sq(5) * Sq(5) == Sq(9, 1)

    assert P(0) * P(5) == P(5)  
    assert P(1) * P(1) == 2*P(2)
    assert P(1) * bP(1) == bP(2) + P(2) * beta
    assert P(2) * P(1) == Adem3.zero()

def test_milnor_mult():
    """
        TODO: fill me in.
    """
    pass
    
def test_adem_to_milnor():
    """
        TODO: fill me in.
    """
    assert type(Sq(2).to_milnor()) == MilnorElement
    assert type(Adem3.P(2).to_milnor()) == MilnorElement
    pass

def test_milnor_to_adem():
    assert type(Milnor2.Sq(2).to_adem()) == AdemElement
    assert type(Milnor3.P(2).to_adem()) == AdemElement

    assert Milnor2.Sq(5).to_adem() == Sq(5)
    assert Milnor2.Sq(0,1).to_adem() == Sq(2,1) + Sq(3)
    assert Milnor2.Sq(1,1).to_adem() == Sq(3,1)
    assert Milnor2.Sq(2,1).to_adem() == Sq(4,1) + Sq(5)
    assert Milnor2.Sq(0,2).to_adem() == Sq(4,2) + Sq(5,1) + Sq(6)
    assert Milnor2.Sq(1,2).to_adem() == Sq(5,2) + Sq(7)
    assert Milnor2.Sq(2,2).to_adem() == Sq(6,2) + Sq(7,1)
    assert Milnor2.Sq(0,0,1).to_adem() == Sq(4,2,1) + Sq(5,2) + Sq(6,1) + Sq(7)
    
    
    assert Milnor3.Q(0).to_adem() == beta
    assert Milnor3.Q(1).to_adem() == P(1) * Adem3.b() + 2 * bP(1)
    assert Milnor3.Q(2).to_adem() == (P(3,1) + 2 * P(4)) * beta + 2 * bP(3,1) + bP(4)
    assert len(Milnor3.Q(3).to_adem()) == 10
    
    assert Milnor3.P(0,1).to_adem() == P(3,1) + 2*P(4)
    assert Milnor3.P(1,1).to_adem() == P(4,1) + P(5)
    assert Milnor3.P(0,2).to_adem() == P(6,2) + 2*P(7,1) + P(8)
    
    # A test case from Sage: "Associativity once failed because of a sign error"
    a,b,c = Milnor3.Q(1), Milnor3.P(3), Milnor3.Q(0,1)
    assert (a*b)*c == a*(b*c)
    
def there_and_back_basis_check():
    for x in [
        Milnor2.Sq(0,1,1), 
        Milnor2.Sq(2,3,1), 
        Milnor3.Q(3)
    ]:
        assert x.to_adem().to_milnor() == x
    
    for x in [
        Sq(8,1),
        Sq(7,3,1), 
        Sq(10,3,1), 
        bP(10,3,1)
    ]:
        assert x.to_milnor().to_adem() == x

def test_adem_antipode():
    antipodes_of_Sqn = [   
        Adem2.unit(),
        Sq(1), 
        Sq(2), 
        Sq(2,1), 
        Sq(3,1) + Sq(4),
        Sq(4,1),
        Sq(4,2),
        Sq(4,2,1),
        Sq(5,2,1) + Sq(6,2) + Sq(7,1) + Sq(8),
        Sq(6,2,1) + Sq(8,1),
        Sq(6,3,1) + Sq(7,3) + Sq(8,2)
    ]
    
    antipodes_of_Pn = [
        Adem3.unit(),
        2*P(1),
        P(2),
        2*P(3),
        P(3,1),
        2*P(4,1) + P(5),
        P(5, 1) + P(6),
        2*P(6, 1),
        P(6, 2),
        2*P(7,2) + P(8,1) + 2*P(9),
    ]
    
    assert beta.antipode() == 2*beta

    
    for n, x in enumerate(antipodes_of_Sqn):
        assert Sq(n).antipode() == x
        
    for n, x in enumerate(antipodes_of_Pn):
        assert P(n).antipode() == x
        
    assert Sq(3).antipode() == Sq(2,1)
    assert Sq(4).antipode() == Sq(3,1) + Sq(4)

    for x in [
        Sq(8,1),
        Sq(7,3,1), 
        Sq(10,3,1), 
        bP(10,3,1)
    ]:
        assert x.antipode().antipode() == x

#
#def implementedByAssignmentLaterInThisFile():
#    assert False, "We implement this by assignment from Vector.linearly_extend_map later in this file Steenrod.py"
#
#
#class AdemElement(Vector):
#    def __init__(self, dictionary, *,  algebra):
#        self.algebra = algebra
#        self.module = algebra
#        super(AdemElement, self).__init__(algebra.p, dictionary)
#    
#    def adem_act(self,v):
#        return self.multiply(v)
#    
#    def multiply(self, v):
#        implementedByAssignmentLaterInThisFile()
#       
#    def to_milnor(self):
#        implementedByAssignmentLaterInThisFile()
#        
#    def basis_degree(self, b):
#        if self.module.generic:
#            result  = 2*(self.p - 1) * sum(b[1::2])
#            result += sum(b[::2])
#        else:
#            result = sum(b)
#        return result
#        
#    def basis_elt_to_string(self, basis_elt):
#        if(self.algebra.generic):
#            result = adem.adem_basis_elt_generic_map(P_fn = lambda P : " P%s" % P, b = " b", basis_elt = basis_elt)
#            return "".join(result)[1:]
#        else:
#            return " ".join(["Sq%s" % s for s in basis_elt]) # This is adem.adem_2_map inlined
#
#AdemElement.multiply = Vector.linearly_extend_map(adem.product, kw_param = "algebra")
#
#class MilnorElement(Vector):
#    def __init__(self, dictionary, *, algebra):
#        self.algebra = algebra
#        self.module = algebra
#        super(MilnorElement, self).__init__(algebra.p, dictionary)
#        
#    def milnor_act(self,v):
#        return self.multiply(v)
#                        
#    def basis_degree(self, b):
#        raise NotImplementedError()
#    
#    def multiply(self, v):        
#        implementedByAssignmentLaterInThisFile()
#            
#    def to_adem(self):
#        implementedByAssignmentLaterInThisFile()
#            
#    def basis_elt_to_string(self, basis_elt):
#        if(self.algebra.generic):
#            Qs = ["Q(%s)" % s for s in basis_elt[0]]
#            Ps = ""
#            if len(basis_elt[1]) > 0:
#                Ps = "P(%s)" % ", ".join([str(s) for s in basis_elt[1]])
#            return " ".join(Qs + [Ps]) or '1'
#        else:
#            Sqs = ""
#            if len(basis_elt) > 0:
#                Sqs = "Sq(%s)" % ", ".join([str(s) for s in basis_elt])
#            return Sqs or '1'
#            
#MilnorElement.multiply = Vector.linearly_extend_map(milnor.product, kw_param = "algebra")
#       
#
#class AdemAlgebra:
#    def __init__(self, p, generic = None):
#        self.p = p
#        if generic is None:
#            generic = p != 2
#        self.generic = generic
#    
#    @staticmethod
#    def getInstance(p, generic = None):
#        if generic is None:
#            generic = p != 2
#        if (p,generic) not in AdemAlgebra.instance_dict:
#            AdemAlgebra.instance_dict[(p,generic)] = AdemAlgebra(p, generic)
#        return AdemAlgebra.instance_dict[(p,generic)]
#    
#    @staticmethod
#    def getInstanceFromAlgebra(algebra):
#        return AdemAlgebra.getInstance(algebra.p, algebra.generic)
#    
#    def getElement(self, d):
#        return AdemElement(d, module = self)
#        
#    def getBasisElement(self, b):
#        return AdemElement({b : 1}, algebra = self)
#    
#    def basis(self, n):
#        return [AdemElement({b:1}, algebra = self) for b in adem.basis(n, algebra = self)]
#    
#    def zero(self):
#        return AdemElement({}, algebra = self)
#
#    def unit(self):
#        return AdemElement({() : 1}, algebra = self)  
#    
#    def b(self):
#        return AdemElement({(1,) : 1}, algebra = self)
#        
#    def b_or_unit(self, epsilon):
#        if self.generic and epsilon:
#            return self.b()
#        else:
#            return self.unit()
#    
#    def Sq(self, n):
#        if self.generic:
#            return AdemElement({(n % 2, n//2, 0) : 1}, algebra = self)
#        else:
#            return AdemElement({(n,) : 1}, algebra = self)
#    
#    def P(self, n):
#        if self.generic:
#            return AdemElement({(0, n, 0) : 1}, algebra = self)
#        else:
#            return AdemElement({(2*n,) : 1}, algebra = self)
#    
#    def bP(self, n):
#        if self.generic:
#            return AdemElement({(1, n, 0) : 1}, algebra = self)
#        else:
#            return AdemElement({(2*n + 1,) : 1}, algebra = self)
#        
#AdemAlgebra.instance_dict = {}
#    
#    
#class MilnorAlgebra(milnor.MinimalMilnorAlgebra):
#    def __init__(self, p, generic = None, profile = None, truncation = None):
#        super(MilnorAlgebra, self).__init__(p, generic, profile, truncation)
#        self.unit_monomial = ((),()) if generic else ()
#    
#    @staticmethod
#    def getInstance(p, generic = None, profile = None, truncation = None):
#        if generic is None:
#            generic = p != 2
#        instance_dict = MilnorAlgebra.instance_dict
#        if (p,generic, profile) not in instance_dict:
#            instance_dict[(p,generic, profile)] = MilnorAlgebra(p, generic, profile, truncation)
#        return instance_dict[(p,generic, profile)]
#    
#    def getInstanceFromAlgebra(algebra):
#        return MilnorAlgebra.getInstance(algebra.p, algebra.generic, algebra.getattr("profile", None), algebra.getattr("truncation", None))
#    
#    def getElement(self, d):
#        return MilnorElement(d, algebra = self)
#
#    def getBasisElement(self, b):
#        return MilnorElement({b : 1}, algebra = self)        
#        
#    def basis(self, n):
#        return [MilnorElement({b:1}, algebra = self) for b in milnor.basis(n, algebra = self)]
#
#    def basis_even(self, n):
#        return [MilnorElement({((), b):1}, algebra = self) for b in milnor.basis_even(n, self.p, self.profile)]
#    
#    def zero(self):
#        return MilnorElement({}, algebra = self)
#    
#    def unit(self):
#        return MilnorElement({ self.unit_monomial : 1 }, algebra = self)
#    
#    def Sq(self, *l):
#        if self.generic:
#            raise NotImplementedError()
#        else:
#            return MilnorElement({ l : 1 }, algebra = self)
#    
#    def Q(self, *l):
#        if self.generic:
#            return MilnorElement({ (l,()) : 1 }, algebra = self)
#        else:
#            return reduce(AdemElement.__mul__, [MilnorElement( { ((0,)* (i) + (1,)) : 1 }, algebra = self) for i in l])
#    
#    def P(self, *l):
#        if self.generic:
#            return MilnorElement({((), l) : 1}, algebra = self)
#        else:
#            return MilnorElement({tuple(2 * i for i in l) : 1}, algebra = self)
#
#MilnorAlgebra.instance_dict = {}
#
#
##@memoized
#def adem_to_milnor_on_basis_2(b):
#    A = MilnorAlgebra.getInstance(2)
#    unit = A.unit()
#    sqs = [A.Sq(i) for i in b]
#    return reduce(lambda a,b : a*b, sqs, unit)
#
##@memoized
#def adem_to_milnor_on_basis_generic(basis_elt, p):
#    A = MilnorAlgebra.getInstance(p)
#    unit = A.unit()
#    sqs = adem.adem_basis_elt_generic_map(P_fn = A.P, b = A.Q(0), basis_elt = basis_elt)
#    return reduce(lambda a,b : a*b, sqs, unit)
#    
#def adem_to_milnor_on_basis(b, *, algebra):
#    if algebra.generic:
#        return adem_to_milnor_on_basis_generic(b, algebra.p)
#    else:
#        return adem_to_milnor_on_basis_2(b)
#
##@memoized
#def milnor_to_adem_on_basis_2(b):
#    t = [0] * len(b)
#    t[-1] = b[-1]
#    for i in range(len(b) - 2, -1, -1):
#        t[i] = b[i] + 2 * t[i + 1]
#    t = tuple(t)
#    x = adem_to_milnor_on_basis_2(t)
#    x.pop(b)
#    result = x.to_adem()
#    result[t] = 1
#    return result
#
##@memoized
#def milnor_to_adem_on_basis_generic(b, *, p):
#    e = b[0]
#    s = b[1]
#    pad_length = max(*e) if e else 0
#    s += (0,) * (pad_length - len(s)) 
#    t = [0,] * (2*len(s) + 1)
#    
#    for i in e:
#        t[2*i] = 1
#            
#    t[-2] = s[-1] + t[-1]
#    for i in range(2, len(s) + 1):
#        t[-2*i] = t[-2*i + 1] + s[-i] + p * t[-2*i + 2]
#    t = tuple(t)
#    x = adem_to_milnor_on_basis_generic(t, p)
#    x.pop(b)
#    x.scale_in_place(-1)
#    result = x.to_adem(output_module = AdemAlgebra.getInstance( p = p, generic = True ))
#    result[t] = 1
#    return result
#
#def milnor_to_adem_on_basis(b, *, algebra):
#    if algebra.generic:
#        return milnor_to_adem_on_basis_generic(b, p = algebra.p)
#    else:
#        return milnor_to_adem_on_basis_2(b)
#
#AdemElement.to_milnor = Vector.linearly_extend_map(adem_to_milnor_on_basis, kw_param = "algebra", 
#                                get_output_module = lambda alg : MilnorAlgebra.getInstance(alg.p, generic = alg.generic))
#MilnorElement.to_adem = Vector.linearly_extend_map(milnor_to_adem_on_basis, kw_param = "algebra",
#                                get_output_module = lambda alg : AdemAlgebra.getInstance(alg.p, generic = alg.generic))
#
#def adem_antipode_on_basis(basis_elt, *, algebra):
#    antipode = algebra.unit()
#    p = algebra.p
#    generic = algebra.generic
#    milnor_alg = MilnorAlgebra.getInstance(p, generic)
#    if not generic:
#        elts = [ Vector.sum(milnor_alg.basis(n)) for n in basis_elt ]
#    else:
#        def P_fn(P):
#            return (-1)**P * Vector.sum( milnor_alg.basis_even( P ) )
#        
#        elts = adem.adem_basis_elt_generic_map(
#            P_fn = P_fn,
#            b = (-1) * milnor_alg.Q(0), 
#            basis_elt = basis_elt
#        )
#    milnor_antipode = reduce( lambda a, b: a*b, elts, milnor_alg.unit() )
#    return milnor_antipode.to_adem()
#
#AdemElement.antipode = Vector.linearly_extend_map( adem_antipode_on_basis, kw_param = "algebra")
#
#
#if __name__ == "__main__":
#    A   = AdemAlgebra(p = 2)
#    A3  = AdemAlgebra(p = 3)
#    Am  = MilnorAlgebra(p = 2)
#    A3m = MilnorAlgebra(p = 3) 
#    
#    print(type(A.Sq(2).to_milnor) == MilnorElement) 
#    print(type(Am.Sq(2).to_adem) == AdemElement)
#    print(type(A3.P(2).to_milnor) == MilnorElement) 
#    print(type(A3m.P(2).to_adem) == AdemElement)
