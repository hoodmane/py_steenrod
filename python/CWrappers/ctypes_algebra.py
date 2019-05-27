from ctypes import *
from ctypes_FpVector import *

#typedef struct Algebra {
#    unsigned long p;
#// Methods:
#    bool (*compute_basis)(struct Algebra* this, unsigned long degree);
#    unsigned long (*get_basis_dimension)(struct Algebra* this, unsigned long degree);
#    int (*multiply_basis_elements)(struct Algebra* this, unsigned long r_degree, unsigned long r, unsigned long s_degree, unsigned long s);
#} Algebra;

class c_Algebra(Structure):
    pass
    
c_Algebra._fields_ = [
        ("p", c_ulong),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong)),
        ("get_basis_dimension", CFUNCTYPE(c_ulong, POINTER(c_Algebra), c_ulong)),
        ("multiply_basis_elements", CFUNCTYPE(c_int, POINTER(c_Algebra), POINTER(c_Vector), c_ulong, c_ulong, c_ulong, c_ulong))
    ]





def wrap_algebra(CSteenrod):
    pass
    
    
    
    
    
    
    
    
    
