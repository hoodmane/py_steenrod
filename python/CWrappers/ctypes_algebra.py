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
        ("algebra", POINTER(c_Algebra)),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong)),
        ("get_basis_dimension", CFUNCTYPE(c_ulong, POINTER(c_Algebra), c_ulong)),
        ("multiply_basis_elements", CFUNCTYPE(c_bool, POINTER(c_Algebra), c_ulong))
    ]

#typedef struct {
#    Module module;
#    unsigned long dimension;
#    unsigned long max_degree;
#    unsigned long * number_of_basis_elements_in_degree;
#    // This goes input_degree --> output_degree --> operation --> input_index --> Vector
#    Vector **** actions;
#} FiniteDimensionalModule;

class c_FiniteDimensionalModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("dimension", c_ulong),
        ("max_degree",c_ulong),
        ("number_of_basis_elements_in_degree", POINTER(c_ulong)),
        ("actions", POINTER(POINTER(POINTER(POINTER(c_Vector)))))
    ]    

#typedef struct {
#    Module module;
#    unsigned long number_of_generators;
#    unsigned long * number_of_generators_in_degree;
#} FreeModule;
class c_FreeModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("number_of_generators", c_ulong),
        ("number_of_generators_in_degree",POINTER(c_ulong)),
    ]    


#typedef struct {
#    FreeModule * source;
#    Module * target;
#    long ** outputs;
#    unsigned long max_kernel_degree;
#    long ** kernel;
#} FreeModuleHomomorphism;
class c_FreeModuleHomomorphism(Structure):
    _fields_ = [
        ("source", POINTER(c_FreeModule)),
        ("target", POINTER(c_Module)),
        ("max_kernel_degree", c_ulong),
        ("kernel", POINTER(POINTER(c_long)))
    ]

#typedef struct {
#    Algebra * algebra;
#    Module * module;
#    FreeModule * resolution_modules;
#    FreeModuleHomomorphism * resolution_differentials;
#} Resolution;

class c_Resolution(Structure):
    _fields_ = [
        ("algebra", c_Algebra),
        ("module", c_Module),
        ("resolution_modules", POINTER(c_FreeModule)),
        ("resolution_differentials", POINTER(c_FreeModule))
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
    
    #FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, unsigned long dimension, unsigned long * generator_degrees);
    CSteenrod.constructFiniteDimensionalModule.argtypes = [POINTER(c_Algebra), c_ulong, POINTER(c_ulong)]
    CSteenrod.constructFiniteDimensionalModule.restype = POINTER(c_FiniteDimensionalModule)
    
    #void freeFiniteDimensionalModule(FiniteDimensionalModule * module);
    CSteenrod.freeFiniteDimensionalModule.argtypes = [POINTER(c_Module)]
    
    #void addActionToFiniteDimensionalModule(FiniteDimensionalModule * module,
    #                                        unsigned long operation_degree, unsigned long operation_idx,
    #                                        unsigned long input_degree, unsigned long input_idx,
    #                                        unsigned long * output
    #);
    CSteenrod.addActionToFiniteDimensionalModule.argtypes = [
        POINTER(c_FiniteDimensionalModule), c_ulong, c_ulong, c_ulong, c_ulong, POINTER(c_ulong)
    ]
    
    
    
    
    
    
    
    
    
