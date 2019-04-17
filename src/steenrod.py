from itertools import product, izip

from memoized import memoized
import adem
import milnor

def add_vectors(v1, v2, p):
    """
        Addition in Fp vector space.
        Add v2 to v1 in place, and reduce resulting keys mod p.
    """
    result = v1 # Don't copy()
    for b in v2:
        if b not in result:
            result[b] = 0
        result[b] += v2.dict[b]
        result[b] = result[b] % p
    return None
    
def linearly_extend_map(p, f):
    """
        V is an Fp vector space with basis B, extend a map f: B1 x ... x Bn -> W to a multilinear map V1 x ... x Vn  -> W.
    """    
    def extension(*args):
        result = {}
        for l in product(*args):
            coeff = 1
            for i in range(len(args)):
                coeff *= args[i][l[i]]
            w = f(*args)
            for b in w:
                w[b] *= coeff
            add_vectors(result, w, p)
        return result
    return extension

@memoized
def get_product(basis, p, generic = None):
    """
        Return the linearly extended product function which takes a pair of 
        vectors in the appropriate basis and returns their product.
        
        basis -- milnor or adem
    """
    if generic is None:
        generic = p != 2
    if generic:
        return linearly_extend_map(p, lambda m1, m2: basis.product_odd(m1, m2, p))
    else:
        return linearly_extend_map(p, lambda m1, m2: type.product_2(m1, m2))

milnor_product_2 = get_product(milnor, 2, False)
adem_product_2   = get_product(adem  , 2, False)

milnor_product_odd = lambda  v, w, p: get_product(milnor, p, True)(v, w)
adem_product_odd   = lambda  v, w, p: get_product(adem,   p, True)(v, w)

@memoized
def adem_to_milnor_on_basis_2(b):
    unit = milnor.unit(2, False)
    mult = get_product(milnor, 2, False)
    sqs = map(milnor.ademSq, b)
    return reduce(mult, sqs, unit)

@memoized
def adem_to_milnor_on_basis_odd(b, p):
    unit = milnor.unit(p, True)
    mult = get_product(milnor, p, True)
    sqs = []
    for idx, j in enumerate(b):
        if idx % 2 == 1:
            if j != 0:
                sqs.append(milnor.Q(0))
        else:
            sqs.append(milnor.ademP(j))
    return reduce(mult, sqs, unit)

def adem_to_milnor_on_basis(b, p, generic = None):
    if generic is None:
        generic = p != 2
    if generic:
        return adem_to_milnor_on_basis_odd(b, p)
    else:
        return adem_to_milnor_on_basis_2(b)

@memoized
def milnor_to_adem_on_basis_2(b):
    t = [0] * len(b)
    t[-1] = b[-1]
    for i in range(len(b) - 2, -1, -1):
        t[i] = b[i] + 2 * t[i + 1]
    x = adem_to_milnor_on_basis_2(tuple(t))
    x.pop(b)
    result = milnor_to_adem_2(x)
    result.set(t,1)
    return result

@memoized
def milnor_to_adem_on_basis_odd(b, p):
    e = b[0]
    s = b[1]
    pad_length = max(*e)
    # Important that we are using tuples here so the input object doesn't
    # get changed!
    s += (0,) * (pad_length - len(s)) 
    t = (0,) * (2*len + 1)
    for i in e:
        t[2*i] = 1
    t[-2] = s[-1] + t[-1]
    for i, idx in izip(range(len(s) - 2, -1, -1), range(len(t) - 2, -1, -2)):
        idx -= 2
        t[idx] = t[idx + 1] + s[i] + p * t[idx + 2]
    x = adem_to_milnor_on_basis_odd(t, p)
    x.pop(b)
    for b1 in x:
        x[b1] *= -1
    result = milnor_to_adem_odd(x, p)
    result.set(t,1)
    return result

def milnor_to_adem_on_basis(b, p, generic):
    """
        Use the upper triangularity of the change of basis matrix from milnor to adem.
        Monks gives us the ordering with respect to which the change of basis is 
        upper triangular and tells us the leading term corresponding to each 
        milnor basis element. 
        
        I guess we wouldn't need to do this if we actually had Fp linear algebra,
        but it's a nice approach.
        
        See Monks paper page 8
    """
    generic = p != 2
    if generic:
        milnor_to_adem_on_basis_odd(b, p)
    else:
        milnor_to_adem_on_basis_2(b)

@memoized
def get_adem_to_milnor(p, generic = None):
    if generic is None:
        generic = p != 2
    if generic:
        return linearly_extend_map(p, lambda b: adem_to_milnor_on_basis_odd(b, p))
    else:
        return linearly_extend_map(p, lambda b: adem_to_milnor_on_basis_2(b))

@memoized
def get_milnor_to_adem(p, generic = None):
    if generic is None:
        generic = p != 2
    if generic:
        return linearly_extend_map(p, lambda b: milnor_to_adem_on_basis_odd(b, p))
    else:
        return linearly_extend_map(p, lambda b: milnor_to_adem_on_basis_2(b))
        
adem_to_milnor_2 = get_adem_to_milnor(2, False)
adem_to_milnor_odd = lambda  v, p: get_adem_to_milnor(p, True)(v)


    
#
#class ModuleElement:
#    __init__(self, p, dict):
#        self.dict = dict
#        self.p = p
#    
#    def __eq__(self, other):
#        return self.dict == other.dict
#    
#    def __add__(self, other):
#        return ModuleElement(self.p, add_vectors(self.p, self.dict.copy(), other.dict))
#            
#class AdemElement(ModuleElement):
#    def __init__(self, p, generic, dict):
#        ModuleElement.__init__(self, p, dict)
#        self.generic = generic
#        
#    def __mult__(self, other):
#        result = {}
#        for (k, l) in product(self.dict, other.dict):
#            coeff = self.dict[k] * other.dict[l]
#            prod = adem.make_mono_admissible(k + l, p, generic)
#            prod = { k : v * coeff for k ,v in prod.items()}
#            add_vectors(self.p, result, prod)
#    
#    def toMilnor(self):
#        pass
#    
#    
#class AdemAlgebra:
#    def __init__(self, p, generic):
#        self.p = p
#        self.generic = generic
#    
#    @staticmethod
#    def getAlgebra(p, generic = None):
#        if generic is None:
#            generic = p != 2
#        if (p,generic) not in instance_dict:
#            AdemAlgebra(self, p, generic)
#        return AdemAlgebra.instance_dict[(p,generic)]
#        
#    def Sq(n):
#    
#    def P(n):
#    
#    def bP(n):
#    
#AdemAlgebra.instance_dict = {}
#        
#
#class MilnorElement(ModuleElement):
#    def __init__(self, p, generic, profile, dict):
#        ModuleElement.__init__(self, p, dict)
#        self.generic = generic
#        self.profile = profile
#    
#    def __mult__(self, other):
#    
#    def __eq__(self, other):
#    
#    def toAdem(self):
#    
#    
#class MilnorAlgebra:
#    def __init__(self, p, generic, profile):
#        if generic is None:
#            generic = p != 2
#        self.p = p
#        self.generic = generic
#        self.profile = profile
#        unit_monomial = ((),()) if generic else ()
#        self.unit = MilnorElement(p, generic, profile, { unit_monomial : 1 })
#    
#    @staticmethod
#    def getAlgebra(p, generic = None, profile = None, truncation = None):
#        if generic is None:
#            generic = p != 2
#        if (p,generic) not in instance_dict:
#            MilnorAlgebra(self, p, generic, profile)
#        return MilnorAlgebra.instance_dict[(p,generic, profile)]
#    
#    def unit():
#        return self.unit
#    
#    def Sq(l):
#    
#    def Q(l):
#    
#    def P(l):
#    
#    
#MilnorAlgebra.instance_dict = {}
