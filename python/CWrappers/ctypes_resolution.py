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
        ("algebra", c_Algebra),
        ("module", c_Module),
        ("resolution_modules", POINTER(c_FreeModule)),
        ("resolution_differentials", POINTER(c_FreeModule))
    ]