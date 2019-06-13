from ctypes import *
from ctypes_algebra import *

# typedef struct {
#     uint degree;
#     uint excess;
#     uint bocksteins;
#     uint P_length;
#     uint *Ps;
# } AdemBasisElement;
class c_AdemBasisElement(Structure):
    _fields_ = [
        ("degree", c_int),
        ("excess", c_int),
        ("bocksteins", c_uint),
        ("P_length", c_uint),
        ("Ps", POINTER(c_uint))
    ]

# typedef struct {
#     uint length;
#     AdemBasisElement **list;
# } AdemBasisElement_list;
class c_AdemBasisElement_list(Structure):
    _fields_ = [
        ("length", c_uint),
        ("list", POINTER(POINTER(c_AdemBasisElement)))
    ]


# typedef struct {
#     Algebra algebra;
#     bool generic;
#     bool unstable;
#     // This will be passed to q_sort and determines our monomial ordering.
#     // return <0 if a should go before b, >0 if a should come after b, and 
#     // 0 if you don't care. It's probably a good idea to impose a total order
#     // for increased consistency.
#     int (*sort_order)(const void *a, const void *b); 
# } AdemAlgebra;
class c_AdemAlgebra(Structure):
    _fields_ = [
        ("algebra", c_Algebra),
        ("generic", c_bool),
        ("unstable", c_bool),
        ("sort_order", CFUNCTYPE(c_int, c_void_p, c_void_p))
    ]


def wrap_adem(CSteenrod):
    # AdemAlgebra *AdemAlgebra_construct(uint p, bool generic, bool unstable);
    CSteenrod.AdemAlgebra_construct.argtypes = [c_uint, c_bool, c_bool]
    CSteenrod.AdemAlgebra_construct.restype = POINTER(c_AdemAlgebra)
    # void AdemAlgebra_free(AdemAlgebra *algebra);s
    CSteenrod.AdemAlgebra_free.argtypes = [c_AdemAlgebra]

    
    # void AdemAlgebra_generateBasis(Algebra *algebra, int max_degree);
    CSteenrod.AdemAlgebra_generateBasis.argtypes = [POINTER(c_Algebra), c_int]
    # void AdemAlgebra_freeBasis(AdemAlgebra *algebra);
    
    # uint AdemAlgebra_getDimension(Algebra *algebra, int degree, int excess);
    CSteenrod.AdemAlgebra_getDimension.argtypes = [POINTER(c_Algebra), c_int, c_int]
    CSteenrod.AdemAlgebra_getDimension.restype = c_uint
    
    # AdemBasisElement_list AdemAlgebra_getBasis(AdemAlgebra *algebra, int degree, int excess);
    CSteenrod.AdemAlgebra_getBasis.argtypes = [POINTER(c_Algebra), c_int, c_int]
    CSteenrod.AdemAlgebra_getBasis.restype = c_AdemBasisElement_list    
    
    # AdemBasisElement AdemAlgebra_basisElement_fromIndex(AdemAlgebra *algebra, int degree, uint idx);
    CSteenrod.AdemAlgebra_basisElement_fromIndex.argtypes = [POINTER(c_AdemAlgebra), c_int, c_uint]
    CSteenrod.AdemAlgebra_basisElement_fromIndex.restype = c_AdemBasisElement
    # uint AdemAlgebra_basisElement_toIndex(AdemAlgebra *public_algebra,  AdemBasisElement *b);
    CSteenrod.AdemAlgebra_basisElement_toIndex.argtypes = [POINTER(c_AdemAlgebra), POINTER(c_AdemBasisElement)]
    CSteenrod.AdemAlgebra_basisElement_toIndex.restype = c_uint   
    
    #void AdemAlgebra_multiply(Algebra * algebra, Vector * result, uint coeff, int r_degree, uint r_index, int s_degree, uint s_index, int excess);
    CSteenrod.AdemAlgebra_multiply.argtypes = [POINTER(c_Algebra), POINTER(c_Vector), c_uint, c_int, c_uint, c_int, c_uint, c_int]

