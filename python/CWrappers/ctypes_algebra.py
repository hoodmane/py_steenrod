from ctypes import *
from ctypes_FpVector import *

# typedef struct {
#     uint length;
#     uint *degrees;
#     uint *indices;
# } FiltrationOneProductList;
class c_FiltrationOneProductList(Structure):
    _fields_ = [
        ("length", c_uint),
        ("degrees", POINTER(c_uint)),
        ("indices", POINTER(c_uint))
    ]

# typedef struct Algebra {
#     uint p;
#     uint max_degree; 
#     char *name;
#     FiltrationOneProductList *product_list; // This determines which indecomposibles have lines drawn for them.
# // Methods:
#     void (*computeBasis)(struct Algebra* this, uint degree);
#     uint (*getDimension)(struct Algebra* this, uint degree, uint excess);
#     void (*multiplyBasisElements)(struct Algebra* this, Vector *result, uint coeff, uint r_degree, uint r, uint s_degree, uint s, uint excess);
# } Algebra;

class c_Algebra(Structure):
    pass
    
c_Algebra._fields_ = [
        ("p", c_uint),
        ("max_degree", c_uint),
        ("name", c_char_p),
        ("product_list", POINTER(c_FiltrationOneProductList)),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong)),
        ("get_basis_dimension", CFUNCTYPE(c_ulong, POINTER(c_Algebra), c_ulong)),
        ("multiply_basis_elements", CFUNCTYPE(c_int, POINTER(c_Algebra), POINTER(c_Vector), c_ulong, c_ulong, c_ulong, c_ulong))
    ]

def wrap_algebra(CSteenrod):
    pass

    
    
    
    
    
    
    
    
    
