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
#     void (*computeBasis)(struct Algebra* this, int degree);
#     uint (*getDimension)(struct Algebra* this, int degree, int excess);
#     void (*multiplyBasisElements)(struct Algebra* this, Vector *result, uint coeff, int r_degree, uint r_idx, int s_degree, uint s_idx, int excess);
# } Algebra;

class c_Algebra(Structure):
    pass
    
c_Algebra._fields_ = [
        ("p", c_uint),
        ("max_degree", c_uint),
        ("name", c_char_p),
        ("product_list", POINTER(c_FiltrationOneProductList)),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Algebra), c_int)),
        ("get_basis_dimension", CFUNCTYPE(c_uint, POINTER(c_Algebra), c_int, c_int)),
        ("multiply_basis_elements", CFUNCTYPE(None, POINTER(c_Algebra), POINTER(c_Vector), c_uint, c_int, c_uint, c_int, c_uint, c_int))
    ]

def wrap_algebra(CSteenrod):
    pass

    
    
    
    
    
    
    
    
    
