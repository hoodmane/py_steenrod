from ctypes import *
from ctypes_algebra import *

# typedef struct Module {
#     uint p;
#     Algebra * algebra;
#     uint type;
#     int min_degree
#     int max_degree; 
# // Methods:
#     bool (*computeBasis)(struct Module* this, int degree);
#     uint (*getDimension)(struct Module* this, int degree);
#     void (*actOnBasis)(struct Module* this, Vector *result, uint coeff, int op_degree,uint op_index, int mod_degree, uint mod_index);
# } Module;

class c_Module(Structure):
    pass

c_Module._fields_ = [
        ("p", c_uint),
        ("algebra", POINTER(c_Algebra)),
        ("type", c_uint),
        ("min_degree", c_int),        
        ("max_degree", c_int),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Module), c_int)),
        ("getBasisDimension", CFUNCTYPE(c_uint, POINTER(c_Module), c_int)),
        ("actOnBasis", CFUNCTYPE(c_uint, 
            POINTER(c_Module), POINTER(c_Vector), 
            c_uint, c_int, c_uint, c_int, c_uint
        ))
    ]

# typedef struct {
#     Module module;
#     int max_basis_degree;
#     uint *graded_dimension;
#     // This goes input_degree --> output_degree --> operation --> input_index --> Vector
#     Vector ****actions;
# } FiniteDimensionalModule;

class c_FiniteDimensionalModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("max_basis_degree", c_int),
        ("graded_dimension", POINTER(c_uint)),
        ("actions", POINTER(POINTER(POINTER(POINTER(c_Vector)))))
    ]    

# typedef struct {
#     Module module;
#     uint *number_of_generators_in_degree;
# } FreeModule;
class c_FreeModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("number_of_generators_in_degree",POINTER(c_uint))
    ]    

# typedef struct {
#     int operation_degree;
#     uint operation_index;
#     int generator_degree;
#     uint generator_index;
# } FreeModuleOperationGeneratorPair;
class c_FreeModuleOperationGeneratorPair(Structure):
    _fields_ = [
        ("operation_degree", c_int),
        ("operation_index" , c_uint),
        ("generator_degree", c_int),
        ("generator_index" , c_uint)
    ]

# typedef struct {
#     uint dimension;
#     uint * column_to_pivot_row;
#     Vector ** kernel;
# } Kernel;

class c_Kernel(Structure):
    _fields_ = [
        ("dimension", c_uint),
        ("column_to_pivot_row", POINTER(c_uint)),
        ("kernel", POINTER(POINTER(c_Vector)))
    ]

# typedef struct {
#     FreeModule *source;
#     Module *target;
#     Vector ***outputs; // degree --> input_idx --> output
#     int max_degree;
#     int max_computed_degree;
#     Matrix **coimage_to_image_isomorphism;
#     Kernel **kernel; // This is redundant with the next module's coimage_to_image_iso
# } FreeModuleHomomorphism;

class c_FreeModuleHomomorphism(Structure):
    _fields_ = [
        ("source", POINTER(c_FreeModule)),
        ("target", POINTER(c_Module)),
        ("outputs", POINTER(POINTER(POINTER(c_Vector)))),
        ("max_degree", c_int),
        ("max_computed_degree", c_int),
        ("coimage_to_image_isomorphism", POINTER(POINTER(c_Matrix))),
        ("kernel", POINTER(POINTER(c_Kernel)))
    ]

def wrap_modules(CSteenrod):
    #FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, int min_degree, int max_degree, uint * generator_degrees);
    CSteenrod.FiniteDimensionalModule_construct.argtypes = [POINTER(c_Algebra), c_int, c_int, POINTER(c_uint)]
    CSteenrod.FiniteDimensionalModule_construct.restype = POINTER(c_FiniteDimensionalModule)
    
    #void freeFiniteDimensionalModule(FiniteDimensionalModule * module);
    CSteenrod.FiniteDimensionalModule_free.argtypes = [POINTER(c_Module)]
    
    #void FiniteDimensionalModule_setAction(FiniteDimensionalModule *module,
    #                                        int operation_degree, uint operation_idx,
    #                                        int input_degree, uint input_idx,
    #                                        Vector *output
    #);
    CSteenrod.FiniteDimensionalModule_setAction.argtypes = [
        POINTER(c_FiniteDimensionalModule), c_int, c_uint, c_int, c_uint, POINTER(c_Vector)
    ]

    # uint FiniteDimensionalModule_getDimension(Module* this, uint degree)
    CSteenrod.FiniteDimensionalModule_getDimension.argtypes = [POINTER(c_FiniteDimensionalModule), c_uint]
    CSteenrod.FiniteDimensionalModule_getDimension.restype  = c_uint

    #void FiniteDimensionalModule_actOnBasis(Module * this, Vector *result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_index)
    CSteenrod.FiniteDimensionalModule_actOnBasis.argtypes = [
        POINTER(c_Module), POINTER(c_Vector), c_uint, c_int, c_uint,c_int, c_uint
    ]


    # FreeModule *FreeModule_construct(Algebra * algebra, int max_degree);
    CSteenrod.FreeModule_construct.argtypes = [POINTER(c_Algebra), c_int]
    CSteenrod.FreeModule_construct.restype  = POINTER(c_FreeModule)

    # void FreeModule_free(FreeModule * module);
    CSteenrod.FreeModule_free.argtypes = [POINTER(c_FreeModule)]

    # void FreeModule_ConstructBlockOffsetTable(FreeModule * M, int degree);
    CSteenrod.FreeModule_ConstructBlockOffsetTable.argtypes = [POINTER(c_FreeModule), c_int]

    # bool FreeModule_computeBasis(Module* this, int degree);
    # uint FreeModule_getDimension(Module* this, int degree);
    CSteenrod.FreeModule_getDimension.argtypes = [POINTER(c_Module), c_int]
    CSteenrod.FreeModule_getDimension.restype = c_uint

    # void FreeModule_actOnBasis(Module * this, Vector * result, uint coeff, int op_degree, uint op_index, int mod_degree, uint mod_idx);
    CSteenrod.FreeModule_actOnBasis.argtypes = [POINTER(c_Module), POINTER(c_Vector), c_uint, c_int, c_uint, c_int, c_uint]

    # uint FreeModule_operationGeneratorToIndex(FreeModule *this, int op_deg, uint op_idx, int gen_deg, uint gen_idx)
    CSteenrod.FreeModule_operationGeneratorToIndex.argtypes \
        = [POINTER(c_FreeModule), c_int, c_uint, c_int, c_uint]
    CSteenrod.FreeModule_operationGeneratorToIndex.restype = c_uint

    # FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, int degree, uint index)
    CSteenrod.FreeModule_indexToOpGen.argtypes = [POINTER(c_FreeModule), c_int, c_uint]
    CSteenrod.FreeModule_indexToOpGen.restype  = c_FreeModuleOperationGeneratorPair

    # FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule * source, Module * target, int max_degree);
    CSteenrod.FreeModuleHomomorphism_construct.argtypes = [POINTER(c_FreeModule), POINTER(c_Module), c_int]
    CSteenrod.FreeModuleHomomorphism_construct.restype = POINTER(c_FreeModuleHomomorphism)
    
    # void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism * f, uint num_gens);
    CSteenrod.FreeModuleHomomorphism_AllocateSpaceForNewGenerators.argtypes = [POINTER(c_FreeModuleHomomorphism), c_uint]

    # void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, int input_degree, uint input_index, Vector *output)
    CSteenrod.FreeModuleHomomorphism_setOutput.argtypes = [POINTER(c_FreeModuleHomomorphism), c_int, c_uint, POINTER(c_Vector)]

    # void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism * f, Vector * result, uint coeff, int input_degree, uint input_index);
    CSteenrod.FreeModuleHomomorphism_applyToBasisElement.argtypes = [POINTER(c_FreeModuleHomomorphism), POINTER(c_Vector), c_uint, c_int, c_uint]

    # void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism * f, Matrix *result, int degree);
    CSteenrod.FreeModuleHomomorphism_getMatrix.argtypes = [POINTER(c_FreeModuleHomomorphism), POINTER(c_Matrix), c_int]

    # Kernel * Kernel_construct(VectorInterface * vectImpl, uint p, uint rows, uint columns);
    CSteenrod.Kernel_construct.argtypes = [c_void_p, c_uint, c_uint, c_uint]
