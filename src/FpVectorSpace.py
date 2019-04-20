import itertools

def implementedByAssignmentLaterInThisFile():
    assert False, "We implement this by assignment from Vector.linearly_extend_map later in this file Steenrod.py"


class Vector(dict):
    def __init__(self, p, d = None):
        self.p = p
        d = d or {}
        super(Vector,self).__init__(d)    
    
    
    def __getitem__(self, idx):
        if idx not in self:
            return 0
        return super(Vector,self).__getitem__(idx)
    
    def __setitem__(self, idx, value):
        value = value % self.p
        if(value == 0):
            self.pop(idx, None)
            return
        super(Vector,self).__setitem__(idx, value % self.p)
    
    def add_in_place(self, v, scale = 1):
        """
            Addition in Fp vector space.
            Add v2 to v1 in place, and reduce resulting keys mod p.
        """
        for b in v:
            self[b] += scale * v[b]
        return None
        
    def scale_in_place(self, c):
        for b in self:
            self[b] *= c
        return None
            
    def __add__(self, v):
        result = self.module.zero()
        result.add_in_place(self)
        result.add_in_place(v)
        return result

    def __sub__(self,v):
        result = Vector(self.p, v)
        result.scale_in_place(-1)
        result.add_in_place(self)
        return result

    def __mul__(self, v):
#        if type(v) == AdemElement:
#            return self.multiply(v)
#        elif callable(getattr(v, "adem_act", None)):
#            return v.adem_act(self)
#        elif type(v) == int:
        if type(v) == int:
            result = self.module.getElementFromDict(self)
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
            result = self.module.getElement(self)
            result.scale_in_place(v)
            return result
        else:
            raise TypeError()

    def new_zero_vector(self):
        return Vector(self.p)        

    @staticmethod
    def sum(args):
        result = args[0].module.zero()
        for v in args:
            result.add_in_place(v)
        return result

    @staticmethod
    def linearly_extend_map(f, kw_param = "module"):
        """
            V is an Fp vector space with basis B, extend a map f: B1 x ... x Bn -> W to a multilinear map V1 x ... x Vn  -> W.
        """    
        def extension(*vecs, output_module = None):
            if output_module is None:
                output_module = vecs[0].module
            if output_module is None:
                raise ValueError()
                
            module_arg = { kw_param : output_module }
                    
            result = output_module.zero()                
                
            for l in itertools.product(*vecs):
                coeff = 1
                for (v, basis_elt) in zip(vecs, l):
                    coeff *= v[basis_elt]
                w = f(*l, **module_arg)
                w = Vector(output_module.p, w)
                result.add_in_place(w, coeff)
            return result
        return extension
        
    @staticmethod
    def tensor_basis_elements(*args):
        return Vector.tensor_symbol.join(args)
    
    def tensor(self, w):
        implementedByAssignmentLaterInThisFile()
        
    def __eq__(self, other):
        raise NotImplementedError()
    
    def basis_degree(self, b):
        """
            Get the degree of a basis vector.
        """
        raise NotImplementedError("basis_degree is overloaded in subclasses, but doesn't have a general definiton.")
    
    def degree(self):
        """
            Returns the degree of an element, or None if the element isn't homogenous.
            Uses basis_degree to compute the degree of basis elements, which is overloaded in subclasses.
        """
        degree_set = set([ self.basis_degree(b) for b in self])
        if len(degree_set) > 1:
            return None
        else:
            return next(iter(degree_set))
    
    def __repr__(self):
        result = []
        for (b,c) in self.items():
            if c == 1:
                coeff = ""
            else:
                coeff = "%s * " % c
            result.append(coeff + self.basis_elt_to_string(b))
        if(len(result) == 0):
            return "0";
        else:
            return "  +  ".join(result)

    def basis_elt_to_string(self, b):
        """
            This is overloaded by subclasses
        """
        return str(b)  

Vector.tensor_symbol = "*"
Vector.tensor = Vector.linearly_extend_map(lambda *args, module = None : { Vector.tensor_basis_elements(*args) : 1 })      
            
