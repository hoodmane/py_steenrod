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


#
#typedef struct {
#    unsigned long degree;
#    unsigned long length;
#    unsigned long *p_part;
#} P_part;

class c_P_part(Structure):
    _fields_ = [
        ("degree", c_ulong),
        ("length", c_ulong),
        ("p_part", POINTER(c_ulong))
    ]

#typedef struct {
#    unsigned long length;
#    P_part *list;
#} P_part_list;

class c_P_part_list(Structure):
    _fields_ = [
        ("length", c_ulong),
        ("list", POINTER(c_P_part))
    ]

#typedef struct {
#    unsigned long degree;
#    unsigned long bit_string;
#} Q_part;

class c_Q_part(Structure):
    _fields_ = [
        ("degree", c_ulong),
        ("bit_string", c_ulong)
    ]

#typedef struct {
#    unsigned long length;
#    Q_part *list;
#} Q_part_list;

class c_Q_part_list(Structure):
    _fields_ = [
        ("length", c_ulong),
        ("list", POINTER(c_Q_part))
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


#
#typedef struct {
#    unsigned long p;
#    bool generic;
#    Profile profile;
#    string name;
#    P_part_list* P_table;
#    P_part_list** P_table_by_P_length;
#    unsigned long P_table_max_degree;
#    Q_part_list* Q_table;
#    unsigned long Q_table_max_tau;
#    MilnorBasisElement_list* basis_table;
#    unsigned long basis_max_degree;
#    khash_t(monomial_index_map) ** basis_name_to_index_map;
#} MilnorAlgebra;

class c_MilnorAlgebra(Structure):
    _fields_ = [
        ("p", c_ulong),
        ("generic", c_bool),
        ("profile", c_Profile),
        ("name", c_char_p),
        ("P_table", POINTER(c_P_part_list)),
        ("P_table_by_P_length",POINTER(POINTER(c_P_part_list))), 
        ("P_table_max_degree", c_ulong),
        ("Q_table", POINTER(c_Q_part_list)),
        ("Q_table_max_tau", c_ulong),
        ("basis_table", POINTER(c_MilnorBasisElement_list)),
        ("basis_max_degree", c_ulong),
        ("basis_name_to_index_map", c_void_p)
    ]

def wrap_milnor_datatypes(CSteenrod):        
    #void milnor_algebra_generate_name(MilnorAlgebra *A);
    CSteenrod.milnor_algebra_generate_name.argtypes = [POINTER(c_MilnorAlgebra)]
         
    # string array_to_string(unsigned long* A, unsigned long length);
    # string milnor_basis_element_to_string(MilnorBasisElement *b);
    # MilnorBasisElement milnor_basis_element_from_string(MilnorAlgebra * algebra, char* elt_string);
    # string milnor_element_to_string(MilnorAlgebra * algebra, Vector * m);
    # string milnor_matrix_to_string(unsigned long** M, unsigned long rows, unsigned long cols);
    
    # void initializeMilnorAlgebraFields(MilnorAlgebra * A);
    CSteenrod.initializeMilnorAlgebraFields.argtypes = [POINTER(c_MilnorAlgebra)]

