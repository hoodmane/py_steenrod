from ctypes import *
from ctypes_modules import *

#typedef struct {
#    Algebra * algebra;
#    Module * module;
#    FreeModule * resolution_modules;
#    FreeModuleHomomorphism * resolution_differentials;
#} Resolution;

class c_Resolution(Structure):
    _fields_ = [
        ("algebra", POINTER(c_Algebra)),
        ("module", POINTER(c_Module)),
        ("max_degree", c_uint),
        ("modules", POINTER(POINTER(c_FreeModule))),
        ("differentials", POINTER(POINTER(c_FreeModuleHomomorphism))),
        ("internal_degree_to_resolution_stage", POINTER(c_uint))
    ]

def wrap_resolution(CSteenrod):
    # Resolution *Resolution_construct(FiniteDimensionalModule *module, uint max_filtration, uint max_degree);
    CSteenrod.Resolution_construct.argtypes = [POINTER(c_FiniteDimensionalModule), c_uint, c_uint]
    CSteenrod.Resolution_construct.restype = POINTER(c_Resolution)
    
    # void Resolution_step(Resolution *resolution, uint homological_degree, uint degree);
    CSteenrod.Resolution_step.argtypes = [POINTER(c_Resolution), c_uint, c_uint]
    
    # void Resolution_resolveThroughDegree(Resolution *res, uint degree);
    CSteenrod.Resolution_resolveThroughDegree.argtypes = [POINTER(c_Resolution), c_uint]