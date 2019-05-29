from ctypes import *
from ctypes_algebra import *

# typedef struct Module {
#     uint p;
#     Algebra * algebra;
# // Methods:
#     bool (*computeBasis)(struct Module* this, uint degree);
#     uint (*getDimension)(struct Module* this, uint degree);
#     void (*actOnBasis)(struct Module* this, Vector *result, uint coeff, uint op_degree,uint op_index, uint mod_degree, uint mod_index);
# } Module;

class c_Module(Structure):
    pass

c_Module._fields_ = [
        ("p", c_uint),
        ("algebra", POINTER(c_Algebra)),
        ("compute_basis",CFUNCTYPE(c_bool, POINTER(c_Module), c_uint)),
        ("getBasisDimension", CFUNCTYPE(c_uint, POINTER(c_Module), c_uint)),
        ("actOnBasis", CFUNCTYPE(c_uint, 
            POINTER(c_Module), POINTER(c_Vector), 
            c_uint, c_uint, c_uint, c_uint, c_uint
        ))
    ]

# typedef struct {
#     Module module;
#     uint dimension;
#     uint max_degree;
#     uint * number_of_basis_elements_in_degree;
#     // This goes input_degree --> output_degree --> operation --> input_index --> Vector
#     Vector **** actions;
# } FiniteDimensionalModule;

class c_FiniteDimensionalModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("dimension", c_uint),
        ("max_degree",c_uint),
        ("number_of_basis_elements_in_degree", POINTER(c_uint)),
        ("actions", POINTER(POINTER(POINTER(POINTER(c_Vector)))))
    ]    

# typedef struct {
#     Module module;
#     uint max_generator_degree;
#     uint max_degree;
#     uint number_of_generators;
#     uint *number_of_generators_in_degree;
# } FreeModule;
class c_FreeModule(Structure):
    _fields_ = [
        ("module", c_Module),
        ("max_generator_degree", c_uint),
        ("max_degree", c_uint),
        ("number_of_generators", c_uint),
        ("number_of_generators_in_degree",POINTER(c_uint))
    ]    

# typedef struct {
#     uint operation_degree;
#     uint operation_index;
#     uint generator_degree;
#     uint generator_index;
# } FreeModuleOperationGeneratorPair;
class c_FreeModuleOperationGeneratorPair(Structure):
    _fields_ = [
        ("operation_degree", c_uint),
        ("operation_index" , c_uint),
        ("generator_degree", c_uint),
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
#     uint max_computed_degree;
#     Matrix **coimage_to_image_isomorphism;
#     Kernel **kernel; // This is redundant with the next module's coimage_to_image_iso
# } FreeModuleHomomorphism;

class c_FreeModuleHomomorphism(Structure):
    _fields_ = [
        ("source", POINTER(c_FreeModule)),
        ("target", POINTER(c_Module)),
        ("outputs", POINTER(POINTER(POINTER(c_Vector)))),
        ("max_computed_degree", c_uint),
        ("coimage_to_image_isomorphism", POINTER(POINTER(c_Matrix))),
        ("kernel", POINTER(POINTER(c_Kernel)))
    ]

def wrap_modules(CSteenrod):
    #FiniteDimensionalModule * constructFiniteDimensionalModule(Algebra * algebra, uint dimension, uint * generator_degrees);
    CSteenrod.FiniteDimensionalModule_construct.argtypes = [POINTER(c_Algebra), c_uint, POINTER(c_uint)]
    CSteenrod.FiniteDimensionalModule_construct.restype = POINTER(c_FiniteDimensionalModule)
    
    #void freeFiniteDimensionalModule(FiniteDimensionalModule * module);
    CSteenrod.FiniteDimensionalModule_free.argtypes = [POINTER(c_Module)]
    
    #void FiniteDimensionalModule_setAction(FiniteDimensionalModule *module,
    #                                        uint operation_degree, uint operation_idx,
    #                                        uint input_degree, uint input_idx,
    #                                        Vector *output
    #);
    CSteenrod.FiniteDimensionalModule_setAction.argtypes = [
        POINTER(c_FiniteDimensionalModule), c_uint, c_uint, c_uint, c_uint, POINTER(c_Vector)
    ]

    # uint FiniteDimensionalModule_getDimension(Module* this, uint degree)
    CSteenrod.FiniteDimensionalModule_getDimension.argtypes = [POINTER(c_FiniteDimensionalModule), c_uint]
    CSteenrod.FiniteDimensionalModule_getDimension.restype  = c_uint

    #void FiniteDimensionalModule_actOnBasis(Module * this, Vector *result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_index)
    CSteenrod.FiniteDimensionalModule_actOnBasis.argtypes = [
        POINTER(c_Module), POINTER(c_Vector), c_uint, c_uint, c_uint,c_uint, c_uint
    ]


    # FreeModule *FreeModule_construct(Algebra * algebra, uint max_generator_degree, uint max_degree);
    CSteenrod.FreeModule_construct.argtypes = [POINTER(c_Algebra), c_uint, c_uint]
    CSteenrod.FreeModule_construct.restype  = POINTER(c_FreeModule)

    # void FreeModule_free(FreeModule * module);
    CSteenrod.FreeModule_free.argtypes = [POINTER(c_FreeModule)]

    # void FreeModule_ConstructBlockOffsetTable(FreeModule * M, uint degree);
    CSteenrod.FreeModule_ConstructBlockOffsetTable.argtypes = [POINTER(c_FreeModule), c_uint]

    # bool FreeModule_computeBasis(Module* this, uint degree);
    # uint FreeModule_getDimension(Module* this, uint degree);
    CSteenrod.FreeModule_getDimension.argtypes = [POINTER(c_Module), c_uint]
    CSteenrod.FreeModule_getDimension.restype = c_uint

    # void FreeModule_actOnBasis(Module * this, Vector * result, uint coeff, uint op_degree, uint op_index, uint mod_degree, uint mod_idx);
    CSteenrod.FreeModule_actOnBasis.argtypes = [POINTER(c_Module), POINTER(c_Vector), c_uint, c_uint, c_uint, c_uint, c_uint]

    # uint FreeModule_operationGeneratorToIndex(FreeModule *this, uint op_deg, uint op_idx, uint gen_deg, uint gen_idx)
    CSteenrod.FreeModule_operationGeneratorToIndex.argtypes \
        = [POINTER(c_FreeModule), c_uint, c_uint, c_uint, c_uint]
    CSteenrod.FreeModule_operationGeneratorToIndex.restype = c_uint

    # FreeModuleOperationGeneratorPair FreeModule_indexToOpGen(FreeModule *this, uint degree, uint index)
    CSteenrod.FreeModule_indexToOpGen.argtypes = [POINTER(c_FreeModule), c_uint, c_uint]
    CSteenrod.FreeModule_indexToOpGen.restype  = c_FreeModuleOperationGeneratorPair

    # FreeModuleHomomorphism *FreeModuleHomomorphism_construct(FreeModule * source, Module * target, uint max_degree);
    CSteenrod.FreeModuleHomomorphism_construct.argtypes = [POINTER(c_FreeModule), POINTER(c_Module), c_uint]
    CSteenrod.FreeModuleHomomorphism_construct.restype = POINTER(c_FreeModuleHomomorphism)
    
    # void FreeModuleHomomorphism_AllocateSpaceForNewGenerators(FreeModuleHomomorphism * f, uint num_gens);
    CSteenrod.FreeModuleHomomorphism_AllocateSpaceForNewGenerators.argtypes = [POINTER(c_FreeModuleHomomorphism), c_uint]

    # void FreeModuleHomomorphism_setOutput(FreeModuleHomomorphism *f, uint input_degree, uint input_index, Vector *output)
    CSteenrod.FreeModuleHomomorphism_setOutput.argtypes = [POINTER(c_FreeModuleHomomorphism), c_uint, c_uint, POINTER(c_Vector)]

    # void FreeModuleHomomorphism_applyToBasisElement(FreeModuleHomomorphism * f, Vector * result, uint coeff, uint input_degree, uint input_index);
    CSteenrod.FreeModuleHomomorphism_applyToBasisElement.argtypes = [POINTER(c_FreeModuleHomomorphism), POINTER(c_Vector), c_uint, c_uint, c_uint]

    # void FreeModuleHomomorphism_getMatrix(FreeModuleHomomorphism * f, Matrix *result, uint degree);
    CSteenrod.FreeModuleHomomorphism_getMatrix.argtypes = [POINTER(c_FreeModuleHomomorphism), POINTER(c_Matrix), c_uint]

    # Kernel * Kernel_construct(VectorInterface * vectImpl, uint p, uint rows, uint columns);
    CSteenrod.Kernel_construct.argtypes = [c_void_p, c_uint, c_uint, c_uint]
