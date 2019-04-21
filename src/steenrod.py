"""
    File: steenrod.py
    Author: Hood Chatham
    Date: 2019/4/21
    License: GPL?

    Defines AdemAlgebra and MilnorAlgebra, which are the Steenrod algebra with 
    the Adem and Milnor basis respecitively. They have elements AdemElement and 
    MilnorElement which extend VectorSpace. 
    
    Mostly this file is just user friendly wrappers for the files adem.py and milnor.py. 
    We define string representations for basis elements and factory methods for 
    producing popular elements. The antipode and basis conversion algorithms are 
    implemented here because they rely on the linearly extended multiplication
    which we don't define until this file.

"""

from functools import reduce
import math

#from memoized import memoized

import adem
import milnor
from FpVectorSpace import Vector, linearextension, linearextension_change_target

class AdemElement(Vector):
    """
        An element of the AdemAlgebra. Defines the Hopf algebra operations on AdemElement.
        These should be produced by the factor class AdemAlgebra. It stores a backreference.
    """
    def __init__(self, dictionary, *, algebra):
        self.algebra = algebra
        self.module = algebra
        self.adem_element = True
        super(AdemElement, self).__init__(algebra.p, dictionary)

    def adem_act(self, v):
        """Compute the product self * v. Implements abstract method of Vector. """
        if not isinstance(v, AdemElement):
            raise TypeError()
        if v.algebra.p != self.algebra.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()
        return v.multiply(self)

    @linearextension
    def multiply(basis_elt_1, basis_elt_2, * ,module):
        """Compute the product self * v """
        return adem.product(basis_elt_1, basis_elt_2, algebra=module)

    @linearextension_change_target(lambda alg: MilnorAlgebra.getInstance(alg.p, generic=alg.generic))
    def to_milnor(basis_elt, *, module):
        """Convert self to a MilnorElement and return that. """
        return milnor_to_adem_on_basis(basis_elt, algebra=module)
    
    @linearextension    
    def antipode(basis_elt, * , module):
        """Compute the antipode of self """
        return adem_antipode_on_basis(basis_elt, algebra=module)
        
    def coproduct(self):
        """Compute the coproduct of self """
        raise NotImplementedError()

    def basis_degree(self, b):
        """Get degree of a basis vector. Implements abstract method of Vector. """
        if self.module.generic:
            result  = 2*(self.p - 1) * sum(b[1::2])
            result += sum(b[::2])
        else:
            result = sum(b)
        return result

    def basis_elt_to_string(self, basis_elt):
        """Get string representation of basis vector. Overrides method in Vector """
        if self.algebra.generic:
            result = adem.adem_basis_elt_generic_map(
                P_fn=lambda P: " P%s" % P,
                b=" b",
                basis_elt=basis_elt
            )
            result = "".join(result)[1:]
        else:
            result = " ".join(["Sq%s" % s for s in basis_elt]) # This is adem.adem_2_map inlined
        return result or "1"

class MilnorElement(Vector):
    """
        An element of the MilnorAlgebra.
        These should be produced by the factor class MilnorAlgebra. It stores a backreference.   
    """
    def __init__(self, dictionary, *, algebra):
        self.algebra = algebra
        self.module = algebra
        self.milnor_element = True
        super(MilnorElement, self).__init__(algebra.p, dictionary)

    def milnor_act(self, v):
        if not isinstance(v, MilnorElement):
            raise TypeError()
        if v.p != self.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()
        return v.multiply(self)

    def basis_degree(self, b):
        raise NotImplementedError()

    @linearextension
    def multiply(basis_element_1, basis_element_2, * , module):
        return milnor.product(basis_element_1, basis_element_2, algebra=module)

    @linearextension_change_target(lambda alg: AdemAlgebra.getInstance(alg.p, generic=alg.generic))
    def to_adem(basis_elt, *, module):
        return milnor_to_adem_on_basis(basis_elt, algebra=module)

    def basis_elt_to_string(self, basis_elt):
        if self.algebra.generic:
            Qs = ["Q(%s)" % s for s in basis_elt[0]]
            Ps = ""
            if basis_elt[1]:
                Ps = "P(%s)" % ", ".join([str(s) for s in basis_elt[1]])
            return " ".join(Qs + [Ps]) or '1'
        else:
            Sqs = ""
            if basis_elt:
                Sqs = "Sq(%s)" % ", ".join([str(s) for s in basis_elt])

            return Sqs or '1'


class AdemAlgebra:
    def __init__(self, p, generic=None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic

    @staticmethod
    def getInstance(p, generic=None):
        if generic is None:
            generic = p != 2
        if (p, generic) not in AdemAlgebra.instance_dict:
            AdemAlgebra.instance_dict[(p, generic)] = AdemAlgebra(p, generic)
        return AdemAlgebra.instance_dict[(p, generic)]

    @staticmethod
    def getInstanceFromAlgebra(algebra):
        return AdemAlgebra.getInstance(algebra.p, algebra.generic)

    def get_element(self, d):
        return AdemElement(d, algebra=self)

    def get_basis_element(self, b):
        return AdemElement({b : 1}, algebra=self)

    def getInadmissiblePairs(self, max_degree):
        p = self.p
        q = 2*(self.p - 1) if self.generic else 1
        P = self.P if self.generic else self.Sq
        for relation_dim in range(2, max_degree // q):
            # We want Pi*b*Pj inadmissible so that means i < p * j + epsilon.
            # relation_dim = i + j so j = relation_dim - i
            # so i < p * (relation_dim - i) + epsilon so i < (p * relation_dim + epsilon) / ( p + 1)
            # We need to round up so that Python includes the last integer
            # if p * relation_dim / (p + 1) is not an integer
            for i in range(1, int(math.ceil((p * relation_dim) / (p+1)))):
                yield (relation_dim, P(i), P(relation_dim - i))
            if self.generic:
                for i in range(1, int(math.ceil((p * relation_dim + 1) / (p+1)))):
                    yield (relation_dim, P(i), self.bP(relation_dim - i))

    def basis(self, n):
        return [AdemElement({b:1}, algebra=self) for b in adem.basis(n, algebra=self)]

    def zero(self):
        return AdemElement({}, algebra=self)

    def unit(self):
        return AdemElement({() : 1}, algebra=self)

    def b(self):
        return AdemElement({(1,) : 1}, algebra=self)

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
        return AdemElement({n : 1}, algebra=self)

    def P(self, *n):
        if not self.generic:
            raise TypeError("Only use A.P if A is generic.")
        # P(0) is the unit
        if len(n) == 1 and n[0] == 0:
            n = []
        l = [0] * (2*len(n) + 1)
        l[1::2] = n
        return AdemElement({tuple(l) : 1}, algebra=self)

    def bP(self, *n):
        if not self.generic:
            raise TypeError("Only use A.bP if A is generic.")
        l = [0] * (2*len(n) + 1)
        l[1::2] = n
        l[0] = 1
        return AdemElement({tuple(l) : 1}, algebra=self)

AdemAlgebra.instance_dict = {}

class MilnorAlgebra(milnor.MinimalMilnorAlgebra):
    def __init__(self, p, generic=None, profile=None, truncation=None):
        super(MilnorAlgebra, self).__init__(p, generic, profile, truncation)
        self.unit_monomial = ((), ()) if generic else ()

    @staticmethod
    def getInstance(p, generic=None, profile=None, truncation=None):
        if generic is None:
            generic = p != 2
        instance_dict = MilnorAlgebra.instance_dict
        if (p, generic, profile) not in instance_dict:
            instance_dict[(p, generic, profile)] = MilnorAlgebra(p, generic, profile, truncation)
        return instance_dict[(p, generic, profile)]

    @staticmethod
    def getInstanceFromAlgebra(algebra):
        return MilnorAlgebra.getInstance(
            algebra.p,
            algebra.generic,
            algebra.getattr("profile", None),
            algebra.getattr("truncation", None)
        )

    def get_element(self, d):
        return MilnorElement(d, algebra=self)

    def get_basis_element(self, b):
        return MilnorElement({b : 1}, algebra=self)

    def basis(self, n):
        return [MilnorElement({b : 1}, algebra=self) for b in milnor.basis(n, algebra=self)]

    def basis_even(self, n):
        return [
            self.get_basis_element(((), b))
            for b in milnor.basis_even(n, self.p, self.profile)
        ]

    def zero(self):
        return MilnorElement({}, algebra=self)

    def unit(self):
        return MilnorElement({self.unit_monomial : 1}, algebra=self)

    def Sq(self, *l):
        if self.generic:
            raise NotImplementedError()
        else:
            return MilnorElement({l : 1}, algebra=self)

    def Q(self, *l):
        if self.generic:
            return MilnorElement({(l, ()) : 1}, algebra=self)
        else:
            return reduce(
                AdemElement.__mul__,
                [self.get_basis_element(((0,)* (i) + (1,))) for i in l]
            )

    def P(self, *l):
        if self.generic:
            return MilnorElement({((), l) : 1}, algebra=self)
        else:
            return MilnorElement({tuple(2 * i for i in l) : 1}, algebra=self)

MilnorAlgebra.instance_dict = {}


#@memoized
def adem_to_milnor_on_basis_2(b):
    A = MilnorAlgebra.getInstance(2)
    unit = A.unit()
    sqs = [A.Sq(i) for i in b]
    return reduce(lambda a, b: a*b, sqs, unit)

#@memoized
def adem_to_milnor_on_basis_generic(basis_elt, p):
    A = MilnorAlgebra.getInstance(p)
    unit = A.unit()
    sqs = adem.adem_basis_elt_generic_map(P_fn=A.P, b=A.Q(0), basis_elt=basis_elt)
    return reduce(lambda a, b: a*b, sqs, unit)

def adem_to_milnor_on_basis(b, *, algebra):
    if algebra.generic:
        return adem_to_milnor_on_basis_generic(b, algebra.p)
    else:
        return adem_to_milnor_on_basis_2(b)

#@memoized
def milnor_to_adem_on_basis_2(b):
    adem_alg = AdemAlgebra.getInstance(p=2, generic=False)
    if not b:
        return adem_alg.unit()
    t = [0] * len(b)
    t[-1] = b[-1]
    for i in range(len(b) - 2, -1, -1):
        t[i] = b[i] + 2 * t[i + 1]
    t = tuple(t)
    x = adem_to_milnor_on_basis_2(t)
    x.pop(b)
    result = x.to_adem(output_module=adem_alg)
    result[t] = 1
    return result

#@memoized
def milnor_to_adem_on_basis_generic(b, *, p):
    (e, s) = b
    adem_alg = AdemAlgebra.getInstance(p=p, generic=True)
    if not e and not s:
        return adem_alg.unit()
    if e == (0,) and not s:
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
    result = x.to_adem(output_module=adem_alg)
    result[t] = 1
    return result

def milnor_to_adem_on_basis(b, *, algebra):
    if algebra.generic:
        return milnor_to_adem_on_basis_generic(b, p=algebra.p)
    else:
        return milnor_to_adem_on_basis_2(b)

def adem_antipode_on_basis(basis_elt, *, algebra):
    p = algebra.p
    generic = algebra.generic
    milnor_alg = MilnorAlgebra.getInstance(p, generic)
    if not generic:
        elts = [Vector.sum(milnor_alg.basis(n)) for n in basis_elt]
    else:
        def P_fn(P):
            return (-1)**P * Vector.sum(milnor_alg.basis_even(P))

        elts = adem.adem_basis_elt_generic_map(
            P_fn=P_fn,
            b=(-1) * milnor_alg.Q(0),
            basis_elt=basis_elt
        )
    # Multiply in the reverse order because antipode is an antihomomorphism
    milnor_antipode = reduce(lambda a, b: b * a, elts, milnor_alg.unit())
    return milnor_antipode.to_adem()


if __name__ == "__main__":
    A   = AdemAlgebra(p=2)
    A3  = AdemAlgebra(p=3)
    Am  = MilnorAlgebra(p=2)
    A3m = MilnorAlgebra(p=3)
