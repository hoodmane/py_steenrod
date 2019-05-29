import json
import itertools
from functools import reduce

from infinity import Infinity
from FpVectorSpace import Vector, linearextension, linearextension_change_target
import steenrod


class FreeModuleElement(Vector):
    def __init__(self, dictionary=None, *,  module):
        self.module = module
        super(FreeModuleElement, self).__init__(module.p, dictionary)

    def adem_act(self, adem_elt):
        if not isinstance(self.module.algebra, steenrod.AdemAlgebra):
            raise TypeError()
        return self.__act_helper(adem_elt)

    def milnor_act(self, milnor_elt):
        if not isinstance(self.module.algebra, steenrod.MilnorAlgebra):
            raise TypeError()
        return self.__act_helper(milnor_elt)

    @linearextension
    def __act_helper(module_basis_elt, alg_basis_elt, *, module):
        result = {}
        (op, gen) = module_basis_elt
        for (k, v) in (module.algebra.get_basis_element(alg_basis_elt) * op).items():
            result[(module.algebra.get_basis_element(k), gen)] = v
        return result


    def basis_degree(self, b):
        return self.module.gens[b[1]] + b[0].degree()
        
    def basis_elt_to_string(self, module, basis_elt):
        result = ""
        if basis_elt[0] != module.algebra.unit():
            result += str(basis_elt[0]) + " "
        result += basis_elt[1]        
        return result
        
    
class FreeModule:
    def __init__(self, *, algebra, name=None):
        self.name = name
        self.algebra = algebra
        self.p = self.algebra.p
        self.gens = {}

    def __repr__(self):
        return "FreeModule(algebra=%s, generators=%s)" % (self.algebra, self.gens)
    
    def basis_degree(self, b):
        return self.gens[b[1]] + b[0].degree()
    
    def add_generator(self, name, degree):
        self.gens[name] = degree
        elt = self.get_generator(name)
        return elt

    def get_generator(self, gen):
        return FreeModuleElement({(self.algebra.unit(), gen) : 1}, module=self)
    
    def get_basis_element(self, op, gen = None):
        if gen:
            opgen = (op, gen)
        else:
            opgen = op
        (op, gen) = opgen
        if gen not in self.gens:
            raise ValueError("%s is not the name of a generator." % gen)
        return FreeModuleElement({(op, gen) : 1}, module=self)

    def get_element(self, d):
        return FreeModuleElement(d, module=self)

    def zero(self):
        return FreeModuleElement({}, module=self)


class ModuleHomomorphism:
    def __init__(self, source, target):
        self.source = source
        self.target = target
        self.algebra = source.algebra
        # self.p = self.algebra.p
        # self.generic = self.algebra.generic
        self.map = {}

    def __repr__(self):
        return "ModuleHomomorphism(source=%s, target=%s)" % (self.source, self.target)

    def add_value(self, input, output):
        if input not in self.source.gens:
            raise ValueError("%s not a generator of source module %s" % (input,self.source))
            
        input_degree = self.source.gens[input]
        output_degree = output.degree()
        if input_degree != output_degree:
            raise ValueError("Target %s not in the appropriate degree. Source %s is in degree %s but target %s is in degree %s" %
                (output, input, input_degree, output, output_degree))
        self.map[input] = output
        
    def __call__(self, input):
        return ModuleHomomorphism.__apply(input, output_module = self.target, map = self.map)
    
    @linearextension
    def __apply(input_basis_vector, *, module, map):
        return input_basis_vector[0] * map[input_basis_vector[1]]


if __name__ == "__main__":
    A = steenrod.MilnorAlgebra(p=2)
    Sq = A.Sq
    M0 = FreeModule(algebra=A)
    x00 = M0.add_generator("x00", 0)
    M1 = FreeModule(algebra=A)
    x11 = M1.add_generator("x11", 1)
    x12 = M1.add_generator("x12", 2)
    x14 = M1.add_generator("x14", 4)
    x18 = M1.add_generator("x18", 8)
    d1 = ModuleHomomorphism(M1, M0)
    d1.add_value("x11", A.Sq(1)*x00)
    d1.add_value("x12", A.Sq(2)*x00)
    d1.add_value("x14", A.Sq(4)*x00)
    d1.add_value("x18", A.Sq(8)*x00)
