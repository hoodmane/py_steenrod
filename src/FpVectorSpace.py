import itertools
import functools

def implementedByAssignmentLaterInThisFile():
    assert False, "We implement this by assignment from Vector.linearly_extend_map later in this file Steenrod.py"

class linearextension_change_target(object):
   """Decorator. Linearly extends a function.
      This version takes an argument which determines the target vector space / module.
      We define the decorator @linearextension to use the vector space that the first argument lives in.
      
      The decorator @linearextension_change_target takes an argument get_output_module.
      The argument get_output_module is a function that takes in the vector space 
      of the first vector and returns the actual target space. I guess the fully 
      general thing would be to get the output vector space as a function of the 
      spaces of all the input vectors but for our applications so far this is sufficient.
   """
   def __init__(self, get_output_module):
      self.get_output_module = get_output_module
      
   def __call__(self, f):
        return extension(self.get_output_module, f)
         
linearextension = linearextension_change_target(None)         
   
class extension(object):
    def __init__(self, get_output_module, func):
        self.get_output_module = get_output_module
        self.func = func
        functools.update_wrapper(self, func)
    
    def __call__(self, *vecs, output_module = None):
        if output_module is None:
            output_module = vecs[0].module
            if self.get_output_module:
                output_module = self.get_output_module(output_module)                    
        if output_module is None:
            raise ValueError()
                
        result = output_module.zero()                
        
        for l in itertools.product(*vecs):
            coeff = 1
            for (v, basis_elt) in zip(vecs, l):
                coeff *= v[basis_elt]
            w = self.func(*l, module=output_module)
            w = Vector(output_module.p, w)
            result.add_in_place(w, coeff)
        return result   
        
    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__
        
    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        if obj is not None:
            return functools.partial(self.__call__, obj)        
        else:
            return self.__call__


class Vector(dict):
    """
        Vector is an abstract class.
        Subclasses should implement the abstract methods:
            basis_degree
            basis_elt_to_string -- default behavior is just to use Str
        
        Many 
        
    """

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
    
    def add_in_place(self, v, coefficient = 1):
        """
            Addition in Fp vector space.
            Add v2 to v1 in place, and reduce resulting keys mod p.
        """
        for b in v:
            self[b] += coefficient * v[b]
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
        if type(v) == int:
            result = self.module.get_element(self)
            result.scale_in_place(v)
            return result
        elif getattr(self, 'adem_element', False) and callable(getattr(v, "adem_act", None)):
            return v.adem_act(self)
        elif getattr(self, 'milnor_element', False) and callable(getattr(v, "milnor_act", None)):
            return v.milnor_act(self)
        else:
            raise TypeError()

    def __rmul__(self, v):
        """
            We could add support here for distinct right and left actions if we want...
        """
        if type(v) == int:
            result = self.module.get_element(self)
            result.scale_in_place(v)
            return result
        elif getattr(self, 'adem_element', False) and callable(getattr(v, "adem_act", None)):
            return v.adem_act(self)
        elif getattr(self, 'milnor_element', False) and callable(getattr(v, "milnor_act", None)):
            return v.milnor_act(self)
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
    def tensor_basis_elements(*args):
        return Vector.tensor_symbol.join(args)
    
    @linearextension
    def tensor(*args, module = None):
        return module.get_basis_element(Vector.tensor_symbol.join(args))
        
#    def __eq__(self, other):
#        raise NotImplementedError()
    
    # abstract method
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

    def basis_elt_to_string(self, basis_elt):
        """
            This is overloaded by subclasses
        """
        return str(basis_elt)  

Vector.tensor_symbol = "*"







