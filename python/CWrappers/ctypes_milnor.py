from ctypes import *
from ctypes_algebra import *

#typedef struct {
#    string name;
#    bool restricted;
#    bool truncated;
#    unsigned long q_part;
#    unsigned long p_part_length;
#    unsigned long* p_part;
#} Profile;
class c_Profile(Structure):
    _fields_ = [
        ("name", c_char_p),
        ("restricted", c_bool),
        ("q_part", c_ulong),
        ("p_part_length", c_ulong),
        ("p_part", POINTER(c_ulong))
    ]

#typedef struct {
#    unsigned long q_degree;
#    unsigned long q_part;
#    unsigned long p_degree;
#    unsigned long p_length;
#    unsigned long *p_part;
#} MilnorBasisElement;
class c_MilnorBasisElement(Structure):
    _fields_ = [
        ("q_degree", c_ulong),
        ("q_part", c_ulong),
        ("p_degree", c_ulong),
        ("p_length", c_ulong),
        ("p_part", POINTER(c_ulong))
    ]

#typedef struct {
#    unsigned long length;
#    MilnorBasisElement *list;
#} MilnorBasisElement_list;
class c_MilnorBasisElement_list(Structure):
    _fields_ = [
        ("length", c_ulong),
        ("list", POINTER(c_MilnorBasisElement))
    ]


#typedef struct {
#    Algebra algebra;
#    unsigned long p;
#    bool generic;
#    Profile profile;
#    string name;
#    unsigned long max_degree;
#} MilnorAlgebra;
class c_MilnorAlgebra(Structure):
    _fields_ = [
        ("p", c_ulong),
        ("generic", c_bool),
        ("profile", c_Profile),
        ("name", c_char_p),
    ]


def wrap_milnor(CSteenrod):
    #// Implemented in milnor_datatypes.c
    #int array_to_string(string buffer, unsigned long* A, unsigned long length);
    #int milnor_basis_element_to_string(string buffer, MilnorBasisElement *b);
    #MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string);
    #
    #// Implemented in milnor_datatypes.c
    #// These methods write a string to a buffer and return the length of the string written.
    #int milnor_element_to_string(string buffer, MilnorAlgebra * algebra, Vector * m);
    #int milnor_matrix_to_string(string buffer, unsigned long** M, unsigned long rows, unsigned long cols);
    #int milnor_basis_element_to_key(string buffer, MilnorBasisElement *b);


    # MilnorAlgebra * constructMilnorAlgebra(unsigned long p, bool generic, Profile *profile);
    CSteenrod.constructMilnorAlgebra.argtypes = [c_ulong, c_bool, POINTER(c_Profile)]
    CSteenrod.constructMilnorAlgebra.restype = POINTER(c_MilnorAlgebra)
    #void freeMilnorAlgebra(MilnorAlgebra *);
    CSteenrod.freeMilnorAlgebra.argtypes = [c_MilnorAlgebra]

    
    #void GenerateMilnorBasis(Algebra * algebra, unsigned long max_degree);
    CSteenrod.GenerateMilnorBasis.argtypes = [POINTER(c_Algebra), c_ulong]
    #void freeMilnorBasis(MilnorAlgebra * algebra);
    
    #unsigned long GetMilnorAlgebraDimension(Algebra * algebra, unsigned long degree);
    CSteenrod.GetMilnorAlgebraDimension.argtypes = [POINTER(c_Algebra), c_ulong]
    CSteenrod.GetMilnorAlgebraDimension.restype = c_ulong
    
    #MilnorBasisElement_list GetMilnorAlgebraBasis(MilnorAlgebra * algebra, unsigned long degree);
    CSteenrod.GetMilnorAlgebraBasis.argtypes = [POINTER(c_Algebra), c_ulong]
    CSteenrod.GetMilnorAlgebraBasis.restype = c_MilnorBasisElement_list    
    
    #MilnorBasisElement GetMilnorBasisElementFromIndex(MilnorAlgebra *algebra, unsigned long degree, unsigned long idx);
    CSteenrod.GetMilnorBasisElementFromIndex.argtypes = [POINTER(c_MilnorAlgebra), c_ulong, c_ulong]
    CSteenrod.GetMilnorBasisElementFromIndex.restype = c_MilnorBasisElement
    #unsigned long GetIndexFromMilnorBasisElement(MilnorAlgebra *algebra,  MilnorBasisElement b);
    CSteenrod.GetIndexFromMilnorBasisElement.argtypes = [POINTER(c_MilnorAlgebra), c_MilnorBasisElement]
    CSteenrod.GetIndexFromMilnorBasisElement.restype = c_ulong   
    
    #void MilnorProduct(Algebra * algebra, Vector * result, unsigned long r_degree, unsigned long r_index, unsigned long s_degree, unsigned long s_index);
    CSteenrod.MilnorProduct.argtypes = [POINTER(c_Algebra), POINTER(c_Vector), c_ulong, c_ulong, c_ulong, c_ulong]

