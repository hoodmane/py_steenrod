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
    
    The only algorithm in here that doesn't come from 

"""

from functools import reduce
import math

from memoized import memoized

import CWrappers
import combinatorics
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
        """Compute the product self * v. Implements abstract method of Vector."""
        if not isinstance(v, AdemElement):
            raise TypeError()
        if v.algebra.p != self.algebra.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()
        return v.multiply(self)

    @linearextension
    def multiply(basis_elt_1, basis_elt_2, *, module):
        """Compute the product self * v."""
        return adem.product(basis_elt_1, basis_elt_2, algebra=module)

    @linearextension_change_target(lambda alg: MilnorAlgebra.getInstanceFromAlgebra(alg))
    def to_milnor(basis_elt, *, module):
        """Convert self to the Milnor basis."""
        return adem_to_milnor_on_basis(basis_elt, p=module.p, generic=module.generic)

    @linearextension
    def antipode(basis_elt, *, module):
        """Compute the antipode of self."""
        return adem_antipode_on_basis(basis_elt, p=module.p, generic=module.generic)

    def coproduct(self):
        """Compute the coproduct of self."""
        raise NotImplementedError()

    def basis_degree(self, b):
        """Get degree of a basis vector. Implements abstract method of Vector."""
        if self.module.generic:
            result = 2*(self.p - 1) * sum(b[1::2])
            result += sum(b[::2])
        else:
            result = sum(b)
        return result

    def basis_elt_to_string(self, module, basis_elt):
        """Get string representation of basis vector. Overrides method in Vector."""
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
    """An element of the MilnorAlgebra.
       These are produced by the factor class MilnorAlgebra and store a
       backreference to the algebra they come from.
    """
    def __init__(self, dictionary, *, algebra):
        self.algebra = algebra
        self.module = algebra
        self.milnor_element = True
        super(MilnorElement, self).__init__(algebra.p, dictionary)

    def milnor_act(self, v):
        """Compute the product self*v. Implements an abstract method of Vector."""
        if not isinstance(v, MilnorElement):
            raise TypeError()
        if v.p != self.p or v.algebra.generic != self.algebra.generic:
            raise TypeError()
        return v.multiply(self)
            

    @linearextension
    def multiply(basis_element_1, basis_element_2, *, module):
        """Compute the product self*v."""
        return milnor.product(basis_element_1, basis_element_2, algebra=module)

    @linearextension_change_target(lambda alg: AdemAlgebra.getInstanceFromAlgebra(alg))
    def to_adem(basis_elt, *, module):
        """Convert self to the Adem basis"""
        return milnor_to_adem_on_basis(basis_elt, p=module.p, generic=module.generic)

    def basis_elt_to_string(self, module, basis_elt):
        """Get string representation of basis vector. Overrides method in Vector."""
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
    """The Steenrod algebra with the Adem basis.
       Primarily this is a factory class for AdemElement, it also keeps track
       of prime and genericness.
    """
    def __init__(self, p, generic=None):
        self.p = p
        if generic is None:
            generic = p != 2
        self.generic = generic
        self.unit_monomial = (0,) if generic else ()

    def __repr__(self):
        result = "AdemAlgebra(p=%s" % self.p
        if self.generic != (self.p != 2):
            result += ", generic=" + self.generic
        result += ")"
        return result
        
    @staticmethod
    def getInstance(p, generic=None):
        """Gets an instance of AdemAlgebra. Same arguments leads to same instance."""
        if generic is None:
            generic = p != 2
        if (p, generic) not in AdemAlgebra.instance_dict:
            AdemAlgebra.instance_dict[(p, generic)] = AdemAlgebra(p, generic)
        return AdemAlgebra.instance_dict[(p, generic)]

    @staticmethod
    def getInstanceFromAlgebra(algebra):
        """Gets an instance of AdemAlgebra with parameters from input algebra (or module)."""
        return AdemAlgebra.getInstance(algebra.p, algebra.generic)

    def get_element(self, d):
        """Turns a dictionary d into an AdemElement."""
        return AdemElement(d, algebra=self)

    def get_basis_element(self, b):
        """Turns a basis element b into an AdemElement."""
        return AdemElement({b : 1}, algebra=self)

    def getInadmissiblePairs(self, max_degree):
        """Returns all triples (a+b, Sqa, Sqb) such that Sqa*Sqb is inadmissible."""
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
        """Returns the list of basis elements in degree n, in no particular order."""
        return [AdemElement({b:1}, algebra=self) for b in adem.basis(n, algebra=self)]

    def zero(self):
        """Returns the zero of the algebra."""
        return AdemElement({}, algebra=self)

    def unit(self):
        """Returns the unit of the algebra."""
        return AdemElement({self.unit_monomial : 1}, algebra=self)

    def b(self):
        """Returns the Bockstein. Works if generic or not."""
        return AdemElement({(1,) : 1}, algebra=self)

    def b_or_unit(self, epsilon):
        """If epsilon is truthy, return the bockstein, else return 1"""
        if self.generic and epsilon:
            return self.b()
        else:
            return self.unit()

    def Sq(self, *n):
        """Returns Sq(i_0) * ... * Sq(i_k). The sequence i need not be admissible."""
        if self.generic:
            raise TypeError("Only use A.Sq if p = 2 and A is not generic.")
        # Sq(0) is the unit
        if len(n) == 1 and n[0] == 0:
            n = ()
        d = adem.make_mono_admissible(n, algebra=self)
        return AdemElement(d, algebra=self)

    def P(self, *n):
        """Returns P(i_0) * ... * P(i_k). The sequence i need not be admissible."""
        if not self.generic:
            raise TypeError("Only use A.P if A is generic.")
        # P(0) is the unit
        if len(n) == 1 and n[0] == 0:
            n = []
        l = [0] * (2*len(n) + 1)
        l[1::2] = n
        d = adem.make_mono_admissible(tuple(l), algebra=self)
        return AdemElement(d, algebra=self)

    def bP(self, *n):
        """Returns b * P(i_0) * ... * P(i_k). The sequence i need not be admissible."""
        if not self.generic:
            raise TypeError("Only use A.bP if A is generic.")
        l = [0] * (2*len(n) + 1)
        l[1::2] = n
        l[0] = 1
        d = adem.make_mono_admissible(tuple(l), algebra=self)
        return AdemElement(d, algebra=self)

AdemAlgebra.instance_dict = {}

class MilnorAlgebra(milnor.MinimalMilnorAlgebra):
    """The Steenrod algebra with the Milnor basis.
       Primarily this is a factory class for MilnorElement, it also keeps track
       of prime, genericness, and profile function.
    """
    def __init__(self, p, generic=None, profile=None, truncation=None):
        super(MilnorAlgebra, self).__init__(p, generic, profile, truncation)
        self.unit_monomial = ((), ()) if generic else ()
    def __repr__(self):
        result = "MilnorAlgebra(p=%s" % (self.p)
        if self.generic != (self.p != 2):
            result += ", generic=%s" % self.generic
        result += ")"
        return result
        

    @staticmethod
    def getInstance(p, generic=None, profile=None, truncation=None):
        """Gets an instance of MilnorAlgebra. Same arguments leads to same instance."""
        if generic is None:
            generic = p != 2
        instance_dict = MilnorAlgebra.instance_dict
        if (p, generic, profile) not in instance_dict:
            instance_dict[(p, generic, profile)] = MilnorAlgebra(p, generic, profile, truncation)
        return instance_dict[(p, generic, profile)]

    @staticmethod
    def getInstanceFromAlgebra(algebra):
        """Gets an instance of MilnorAlgebra with same properties as input algebra."""
        return MilnorAlgebra.getInstance(
            algebra.p,
            algebra.generic,
            getattr(algebra,"profile", None),
            getattr(algebra,"truncation", None)
        )

    def basis_q_degree(self,b):
        if not self.generic:
            return 0
        tau_degrees = combinatorics.tau_degrees(10000, p = self.p);
        Qs = b[0]
        result = 0            
        for i in Qs:
            result += tau_degrees[i]                
        return result    
        
    def basis_p_degree(self, b):
        if self.generic:    
            Ps = b[1]
            q = 2 * (self.p - 1)
        else:
            Ps = b
            q = 1
        xi_degrees = combinatorics.xi_degrees(10000, p = self.p);            
        result = 0               
        for i, exponent in enumerate(Ps):
            result += xi_degrees[i] * exponent * q
        return result    

    def basis_degree(self, b):
        """Get degree of a basis vector. Implements abstract method of Vector."""
        xi_degrees = combinatorics.xi_degrees(10000, p = self.p);
        if self.generic:
            tau_degrees = combinatorics.tau_degrees(10000, p = self.p);
            Qs = b[0]
            Ps = b[1]
            result = 0
            for i in Qs:
                result += tau_degrees[i]
            for i, exponent in enumerate(Ps):
                result += xi_degrees[i] * exponent * 2 * (self.p - 1)                
            return result
        else:
            result = 0
            for i, exponent in enumerate(b):
                result += xi_degrees[i] * exponent
            return result

    def get_element(self, d):
        """Makes a MilnorElement from a dictionary.
           No checking to determine if this is a valid element...
        """
        return MilnorElement(d, algebra=self)

    def get_basis_element(self, b):
        """Makes a MilnorElement from a basis vector.
           No checking to determine if this is a valid basis element...
        """
        return MilnorElement({b : 1}, algebra=self)

    def basis(self, n):
        """Returns the Milnor basis in degree n"""
        return [MilnorElement({b : 1}, algebra=self) for b in milnor.basis(n, algebra=self)]

    def basis_even(self, n):
        """Returns the even part of the Milnor basis in degree n"""
        return [
            self.get_basis_element(((), b))
            for b in milnor.basis_even(n, self.p, self.profile)
        ]

    def zero(self):
        """Returns the zero of the algebra"""
        return MilnorElement({}, algebra=self)

    def unit(self):
        """Returns the unit of the algebra"""
        return MilnorElement({self.unit_monomial : 1}, algebra=self)

    def Sq(self, *l):
        """Returns the Milnor basis element Sq(i_1, ..., i_k)
           By definition this is the element dual to Prod xi_j ^(i_j)
        """
        if self.generic:
            raise NotImplementedError()
        else:
            return MilnorElement({l : 1}, algebra=self)

    def Q(self, *l):
        """Returns the milnor element Q(i_1) ... Q(i_n).
           By definition, this is the element dual to tau_{i_1} ... tau_{i_k}.
           The i's should be distinct and increasing (this might be relaxed later).
        """
        if self.generic:
            return MilnorElement({(l, ()) : 1}, algebra=self)
        else:
            return reduce(
                AdemElement.__mul__,
                [self.get_basis_element(((0,)* (i) + (1,))) for i in l]
            )

    def P(self, *l):
        """Returns the milnor basis element P(i_1, ..., i_n).
           By definition this is the element dual to Prod xi_j ^(i_j)
        """
        if self.generic:
            return MilnorElement({((), l) : 1}, algebra=self)
        else:
            return MilnorElement({tuple(2 * i for i in l) : 1}, algebra=self)


MilnorAlgebra.instance_dict = {}

@memoized
def adem_to_milnor_on_basis(basis_elt, *, p, generic):
    """Convert an Adem basis element to the Milnor basis.
       This is easy because Sqi is a Milnor basis element and Adem basis elements
       are products of Sqi's so we just move over the Sqi's and then multiply
       in the MilnorAlgebra.
    """
    A = MilnorAlgebra.getInstance(p=p, generic=generic)
    unit = A.unit()
    if generic:
        sqs = adem.adem_basis_elt_generic_map(P_fn=A.P, b=A.Q(0), basis_elt=basis_elt)
    else:
        sqs = [A.Sq(i) for i in basis_elt]
    return reduce(lambda a, b: a*b, sqs, unit)

@memoized
def milnor_to_adem_on_basis_2(b):
    """Convert a Milnor basis element to the Adem basis when p=2.
       See milnor_to_adem_on_basis for some info about the algorithm.
    """
    adem_alg = AdemAlgebra.getInstance(p=2, generic=False)
    if not b:
        return adem_alg.unit()
    t = [0] * len(b)
    t[-1] = b[-1]
    for i in range(len(b) - 2, -1, -1):
        t[i] = b[i] + 2 * t[i + 1]
    t = tuple(t)
    x = adem_to_milnor_on_basis(t, p=2, generic=False)
    x.pop(b)
    result = x.to_adem(output_module=adem_alg)
    result[t] = 1
    return result

@memoized
def milnor_to_adem_on_basis_generic(b, *, p):
    """Convert a Milnor basis element to the Adem basis when p=2.
       See milnor_to_adem_on_basis for some info about the algorithm.
    """
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
    x = adem_to_milnor_on_basis(t, p=p, generic=True)
    x.pop(b)
    x.scale_in_place(-1)
    result = x.to_adem(output_module=adem_alg)
    result[t] = 1
    return result

def milnor_to_adem_on_basis(b, *, p, generic):
    """This is the algorithm for inverting a filtered isomorphism (AKA triangular matrix)
       given an ordering on the basis and a function from one basis to the other.

       The computation of the ordering / function is determined by Monks:
       Change of Basis, Monomial Relations, and Pst Bases for the Steenrod Algebra
       https://monks.scranton.edu/files/pubs/bases.pdf
       The formula we use here appears on page 8.

       Monk's paper only treats the case p = 2, I worked out the straightforward
       generalization to the generic case.
    """
    if generic:
        return milnor_to_adem_on_basis_generic(b, p=p)
    else:
        return milnor_to_adem_on_basis_2(b)

def adem_antipode_on_basis(basis_elt, *, p, generic):
    """Compute the antipode on an Adem basis element.

       Stolen from Sage with some neatening. Here's the Sage documentation:
       
       ALGORITHM: according to Milnor, the antipode of Sq(n) is the sum of all the
       Milnor basis elements in dimension n. We use this together with the fact that
       the antipode is an antihomomorphism.

       At odd primes, the antipode of P(n) is the sum of the Milnor P basis
       elements in dimension n*2(p-1), multiplied by `(-1)^n`, and the antipode
       of \beta = Q_0 is -Q_0.
    """
    milnor_alg = MilnorAlgebra.getInstance(p, generic)
    if not generic:
        elts = [Vector.sum(milnor_alg.basis(n)) for n in basis_elt]
    else:
        def P_fn(n):
            return (-1)**n * Vector.sum(milnor_alg.basis_even(n))

        elts = adem.adem_basis_elt_generic_map(
            P_fn=P_fn,
            b=(-1) * milnor_alg.Q(0),
            basis_elt=basis_elt
        )
    # Multiply in the reverse order because antipode is an antihomomorphism
    milnor_antipode = reduce(lambda a, b: b * a, elts, milnor_alg.unit())
    return milnor_antipode.to_adem()


def nilpotence_order(op):
    x = op
    i = 1
    while x:
        x *= op
        i += 1
    return i

if __name__ == "__main__":
    A = AdemAlgebra(p=2)
    A3 = AdemAlgebra(p=3)
    Am = MilnorAlgebra(p=2)
    A3m = MilnorAlgebra(p=3)
    
    P2 = Am.P(0,2)
    print(nilpotence_order(P2))
    
