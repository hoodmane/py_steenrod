from ctypes import *
from ctypes_modules import *

# typedef struct Resolution_s {
#     Algebra *algebra;
#     Module *module;
#     void (*addClass)(uint hom_deg, uint int_deg, char *cocycle_name);
#     void (*addStructline)(
#         uint source_hom_deg, uint source_int_deg, uint source_idx, 
#         uint target_hom_deg, uint target_int_deg, uint target_idx
#     );
#     uint max_degree;
#     FreeModule **modules; // The index into resolution_modules is homological_degree + 1.
#     FreeModuleHomomorphism **differentials;// Each differential has source the module with the same index in resolution_modules
#     int *internal_degree_to_resolution_stage;       // Records how far we've resolved in each degree (homological_degree + 1)
# } Resolution;

class c_Resolution(Structure):
    _fields_ = [
        ("algebra", POINTER(c_Algebra)),
        ("module", POINTER(c_Module)),
        ("addClass", c_void_p),
        ("addStructline", c_void_p),
        ("max_degree", c_uint),
        ("modules", POINTER(POINTER(c_FreeModule))),
        ("differentials", POINTER(POINTER(c_FreeModuleHomomorphism))),
        ("internal_degree_to_resolution_stage", POINTER(c_uint))
    ]

def wrap_resolution(CSteenrod):
    # Resolution *Resolution_construct(FiniteDimensionalModule *module, uint max_filtration, uint max_degree);
    CSteenrod.Resolution_construct.argtypes = [POINTER(c_FiniteDimensionalModule), c_uint, c_uint]
    CSteenrod.Resolution_construct.restype = POINTER(c_Resolution)
    
    CSteenrod.testResolution.restype = POINTER(c_Resolution)

    # void Resolution_step(Resolution *resolution, uint homological_degree, uint degree);
    CSteenrod.Resolution_step.argtypes = [POINTER(c_Resolution), c_uint, c_uint]
    
    # void Resolution_resolveThroughDegree(Resolution *res, uint degree);
    CSteenrod.Resolution_resolveThroughDegree.argtypes = [POINTER(c_Resolution), c_uint]