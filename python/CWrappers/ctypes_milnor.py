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
        ("name", c_char_p),
        ("restricted", c_bool),
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
        ("q_degree", c_uint),
        ("q_part", c_uint),
        ("p_degree", c_uint),
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


#typedef struct {
#    Algebra algebra;
#    uint p;
#    bool generic;
#    Profile profile;
#    string name;
#    uint max_degree;
#} MilnorAlgebra;
class c_MilnorAlgebra(Structure):
    _fields_ = [
        ("p", c_uint),
        ("generic", c_bool),
        ("profile", c_Profile),
        ("name", c_char_p),
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


    # MilnorAlgebra * constructMilnorAlgebra(uint p, bool generic, Profile *profile);
    CSteenrod.constructMilnorAlgebra.argtypes = [c_uint, c_bool, POINTER(c_Profile)]
    CSteenrod.constructMilnorAlgebra.restype = POINTER(c_MilnorAlgebra)
    #void freeMilnorAlgebra(MilnorAlgebra *);
    CSteenrod.freeMilnorAlgebra.argtypes = [c_MilnorAlgebra]

    
    #void GenerateMilnorBasis(Algebra * algebra, uint max_degree);
    CSteenrod.GenerateMilnorBasis.argtypes = [POINTER(c_Algebra), c_uint]
    #void freeMilnorBasis(MilnorAlgebra * algebra);
    
    #uint GetMilnorAlgebraDimension(Algebra * algebra, uint degree);
    CSteenrod.GetMilnorAlgebraDimension.argtypes = [POINTER(c_Algebra), c_uint]
    CSteenrod.GetMilnorAlgebraDimension.restype = c_uint
    
    #MilnorBasisElement_list GetMilnorAlgebraBasis(MilnorAlgebra * algebra, uint degree);
    CSteenrod.GetMilnorAlgebraBasis.argtypes = [POINTER(c_Algebra), c_uint]
    CSteenrod.GetMilnorAlgebraBasis.restype = c_MilnorBasisElement_list    
    
    #MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, uint degree, uint idx);
    CSteenrod.GetMilnorBasisElementFromIndex.argtypes = [POINTER(c_MilnorAlgebra), c_uint, c_uint]
    CSteenrod.GetMilnorBasisElementFromIndex.restype = c_MilnorBasisElement
    #uint GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);
    CSteenrod.GetIndexFromMilnorBasisElement.argtypes = [POINTER(c_MilnorAlgebra), c_MilnorBasisElement]
    CSteenrod.GetIndexFromMilnorBasisElement.restype = c_uint   
    
    #void MilnorProduct(Algebra * algebra, Vector * result, uint coeff, uint r_degree, uint r_index, uint s_degree, uint s_index);
    CSteenrod.MilnorProduct.argtypes = [POINTER(c_Algebra), POINTER(c_Vector), c_uint, c_uint, c_uint, c_uint, c_uint]

