from ctypes import *
from ctypes_algebra import *

#typedef struct {
#    string name;
#    bool restricted;
#    bool truncated;
#    uint q_part;
#    uint p_part_length;
#    uint* p_part;
#} Profile;
class c_Profile(Structure):
    _fields_ = [
        ("generic", c_bool),
        ("restricted", c_bool),
        ("truncated", c_bool),
        ("q_part", c_uint),
        ("p_part_length", c_uint),
        ("p_part", POINTER(c_uint))
    ]

#typedef struct {
#    uint q_degree;
#    uint q_part;
#    uint p_degree;
#    uint p_length;
#    uint *p_part;
#} MilnorBasisElement;
class c_MilnorBasisElement(Structure):
    _fields_ = [
        ("q_degree", c_int),
        ("q_part", c_uint),
        ("p_degree", c_int),
        ("p_length", c_uint),
        ("p_part", POINTER(c_uint))
    ]

#typedef struct {
#    uint length;
#    MilnorBasisElement *list;
#} MilnorBasisElement_list;
class c_MilnorBasisElement_list(Structure):
    _fields_ = [
        ("length", c_uint),
        ("list", POINTER(c_MilnorBasisElement))
    ]


# typedef struct {
#     Algebra algebra;
#     bool generic;
#     Profile profile;
#     int (*sort_order)(const void *a, const void *b); 
# } MilnorAlgebra;
class c_MilnorAlgebra(Structure):
    _fields_ = [
        ("algebra", c_Algebra),
        ("generic", c_bool),
        ("profile", c_Profile),
        ("sort_order", CFUNCTYPE(c_int, c_void_p, c_void_p))
    ]


def wrap_milnor(CSteenrod):
    #// Implemented in milnor_datatypes.c
    #int array_to_string(string buffer, uint* A, uint length);
    #int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b);
    #MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string);
    #
    #// Implemented in milnor_datatypes.c
    #// These methods write a string to a buffer and return the length of the string written.
    #int milnor_element_to_string(string buffer, MilnorAlgebra * algebra, Vector * m);
    #int milnor_matrix_to_string(string buffer, uint** M, uint rows, uint cols);
    #int milnor_basis_element_to_key(string buffer, MilnorBasisElement *b);


    # MilnorAlgebra * MilnorAlgebra_construct(uint p, bool generic, Profile *profile);
    CSteenrod.MilnorAlgebra_construct.argtypes = [c_uint, c_bool, POINTER(c_Profile)]
    CSteenrod.MilnorAlgebra_construct.restype = POINTER(c_MilnorAlgebra)
    #void MilnorAlgebra_free(MilnorAlgebra *);
    CSteenrod.MilnorAlgebra_free.argtypes = [c_MilnorAlgebra]

    
    #void MilnorAlgebra_generateBasis(Algebra * algebra, int max_degree);
    CSteenrod.MilnorAlgebra_generateBasis.argtypes = [POINTER(c_Algebra), c_int]
    #void freeMilnorBasis(MilnorAlgebra * algebra);
    
    #uint MilnorAlgebra_getDimension(Algebra * algebra, int degree, int excess);
    CSteenrod.MilnorAlgebra_getDimension.argtypes = [POINTER(c_Algebra), c_int, c_int]
    CSteenrod.MilnorAlgebra_getDimension.restype = c_uint
    
    #MilnorBasisElement_list MilnorAlgebra_getBasis(MilnorAlgebra * algebra, int degree, int excess);
    CSteenrod.MilnorAlgebra_getBasis.argtypes = [POINTER(c_Algebra), c_int, c_int]
    CSteenrod.MilnorAlgebra_getBasis.restype = c_MilnorBasisElement_list    
    
    #MilnorBasisElement MilnorAlgebra_basisElement_fromIndex(MilnorAlgebra *algebra, int degree, uint idx);
    CSteenrod.MilnorAlgebra_basisElement_fromIndex.argtypes = [POINTER(c_MilnorAlgebra), c_int, c_uint]
    CSteenrod.MilnorAlgebra_basisElement_fromIndex.restype = c_MilnorBasisElement
    #uint MilnorAlgebra_basisElement_toIndex(MilnorAlgebra *algebra,  MilnorBasisElement b);
    CSteenrod.MilnorAlgebra_basisElement_toIndex.argtypes = [POINTER(c_MilnorAlgebra), c_MilnorBasisElement]
    CSteenrod.MilnorAlgebra_basisElement_toIndex.restype = c_uint   
    
    #void MilnorAlgebra_multiply(Algebra * algebra, Vector * result, uint coeff, int r_degree, uint r_index, int s_degree, uint s_index, int excess);
    CSteenrod.MilnorAlgebra_multiply.argtypes = [POINTER(c_Algebra), POINTER(c_Vector), c_uint, c_int, c_uint, c_int, c_uint, c_int]

