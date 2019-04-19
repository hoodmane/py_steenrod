import itertools

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
        super(Vector,self).__setitem__(idx, value % self.p)
    
    def add_in_place(self, v):
        """
            Addition in Fp vector space.
            Add v2 to v1 in place, and reduce resulting keys mod p.
        """
        for b in v:
            self[b] += v[b]
        return None
        
    def scale_in_place(self, c):
        for b in self:
            self[b] *= c
        return None
            
    def __add__(self, v):
        result = Vector(self.p, self)
        result.add_vector(v)
        return result

    def __sub__(self,v):
        result = Vector(self.p, v)
        result.scale_in_place(-1)
        result.add_in_place(self)
        return result

    def new_zero_vector(self):
        return Vector(self.p)        

    @staticmethod
    def sum(args):
        result = Vector(args[0].p)
        for v in args:
            result.add_vector(v)
        return result

    @staticmethod
    def linearly_extend_map(f):
        """
            V is an Fp vector space with basis B, extend a map f: B1 x ... x Bn -> W to a multilinear map V1 x ... x Vn  -> W.
        """    
        def extension(*args):
            alg = args[0].algebra
            result = args[0].algebra.zero()
            for l in itertools.product(*args):
                coeff = 1
                for i in range(len(args)):
                    coeff *= args[i][l[i]]
                w = f(*l, algebra = alg)
                w = Vector(alg.p, w)
                w.scale_in_place(coeff)
                
                result.add_in_place(w)
            return result
        return extension
        
    def __eq__(self, other):
        raise NotImplementedError()
    
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
        return str(b)        
            
