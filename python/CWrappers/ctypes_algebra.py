from ctypes import *

#typedef struct {
#    unsigned long p;
#    unsigned long degree;
#    unsigned long dimension;
#    long* vector;
#} MilnorElement;

class c_Vector(Structure):
    _fields_ = [
        ("p", c_ulong),
        ("degree", c_ulong),
        ("dimension", c_ulong),
        ("vector", POINTER(c_long))
    ]
    
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


#typedef struct Module {
#    unsigned long p;
#    Algebra * algebra;
#// Methods:
#    bool (*compute_basis)(struct Module* this, unsigned long degree);
#    unsigned long (*get_basis_dimension)(struct Module* this, unsigned long degree);
#    int (*act_on_basis)(struct Module* this, unsigned long op_degree, unsigned long op_index, unsigned long mod_degree, unsigned long mod_index);
#} Module;

class c_Module(Structure):
    _fields_ = [
        ("p", c_ulong),
        ("algebra", c_Algebra),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong)),
        ("get_basis_dimension", CFUNCTYPE(c_ulong, POINTER(c_Algebra), c_ulong)),
        ("multiply_basis_elements", CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong))
    ]

def wrap_algebra(CSteenrod):        
    # Vector * allocateVector(unsigned long p, unsigned long degree, unsigned long dimension);
    CSteenrod.allocateVector.argtypes = [c_ulong, c_ulong, c_ulong]
    CSteenrod.allocateVector.restype = POINTER(c_Vector)
    # void freeVector(Vector * elt);
    CSteenrod.freeVector.argtypes = [POINTER(c_Vector)]
         
    # void addBasisElementToVector(Vector * elt, unsigned long idx, long coeff);
    # void addVector(Vector * target, Vector * source);
    # void scaleVector(Vector *, long);
    # void assignVector(Vector * target, Vector * source);
