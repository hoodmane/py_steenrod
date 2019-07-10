/**
 * FpVector.h
 */

// void initializePrime(uint p);
let cinitializePrime = cwrap("initializePrime", 'void', ['number']);

// Vector *Vector_construct(uint p, uint dimension, uint offset);
let cVector_construct = cwrap("Vector_construct", 'pointer', ['number', 'number', 'number']);

// void Vector_free(Vector *v); 
let cVector_free = cwrap("Vector_free", 'void', ['pointer']);

// void Vector_pack(Vector *target, uint *source);
let cVector_pack = cwrap("Vector_pack", 'void', ['pointer', 'pointer']);

// void Vector_unpack(uint * target, Vector * source);
let cVector_unpack = cwrap("Vector_unpack", 'void', ['pointer', 'pointer']);

let cVector_print = cwrap("Vector_print", 'void', ['pointer', 'pointer']);

/**
 * algebra.h
 */
// void algebra_computeBasis_function(Algebra *algebra, int degree);
let cAlgebra_computeBasis = cwrap('Algebra_computeBasis_function', 'void', ['pointer', 'number']);

// uint algebra_getDimension_function(Algebra *algebra, int degree, int excess);
let cAlgebra_getDimension = cwrap('Algebra_getDimension_function', 'number', ['pointer', 'number', 'number']);

// void algebra_multiplyBasisElements_function(Algebra *algebra, Vector *result, uint coeff, int r_deg, uint r_idx, int s_deg, uint s_idx, int excess);
let cAlgebra_multiplyBasisElements = cwrap('Algebra_multiplyBasisElements_function', 'void', ['pointer', 'pointer', 'number', 'number', 'number', 'number', 'number', 'number']);

// uint algebra_basisElementToString_function(Algebra *algebra, char *result, int degree, uint idx);
let cAlgebra_basisElementToString_function = cwrap('Algebra_basisElementToString_function', 'number', ['pointer', 'pointer', 'number', 'number']);


/**
 * adem.h
 */


// AdemAlgebra *AdemAlgebra_construct(uint p, bool generic, bool unstable);
let cAdemAlgebra_construct = Module.cwrap("AdemAlgebra_construct", 'pointer', ['number', 'bool', 'bool']);

// void AdemAlgebra_free(AdemAlgebra *algebra);
let cAdemAlgebra_free = Module.cwrap("AdemAlgebra_free", 'void', ['pointer']);

// void AdemAlgebra_generateBasis(Algebra *this, int max_degree);
let cAdemAlgebra_generateBasis = Module.cwrap("AdemAlgebra_generateBasis", 'bool', ['pointer', 'number']);

// uint AdemAlgebra_getDimension(Algebra *this, int degree, int excess);
let cAdemAlgebra_getDimension = Module.cwrap("AdemAlgebra_getDimension", 'number', ['pointer', 'number', 'number']);

// uint AdemAlgebra_element_toString(char *buffer, AdemAlgebra *algebra, int degree, Vector *m);
let cAdemAlgebra_element_toString = Module.cwrap("AdemAlgebra_element_toString", 'number', ['pointer', 'pointer', 'number', 'pointer']);

// AdemBasisElement *AdemAlgebra_basisElement_fromIndex(AdemAlgebra *public_algebra, int degree, uint index);
let cAdemAlgebra_basisElement_fromIndex = Module.cwrap("AdemAlgebra_basisElement_fromIndex", 'pointer', ['pointer', 'number', 'number']);

// AdemBasisElement *AdemAlgebra_basisElement_construct(uint degree, uint P_length, uint *Ps, uint bocksteins);
let cAdemAlgebra_basisElement_construct = Module.cwrap("AdemAlgebra_basisElement_construct", 'pointer', ['number','number', 'pointer', 'number']);

// uint AdemAlgebra_basisElement_getPlength(AdemBasisElement *b);
let cAdemAlgebra_basisElement_getPlength = Module.cwrap("AdemAlgebra_basisElement_getPlength", 'number', ['pointer']);

// uint *AdemAlgebra_basisElement_getPs(AdemBasisElement *b);
let cAdemAlgebra_basisElement_getPs = Module.cwrap("AdemAlgebra_basisElement_getPs", 'pointer', ['pointer']);

// uint AdemAlgebra_basisElement_getBocksteins(AdemBasisElement *b);\
let cAdemAlgebra_basisElement_getBocksteins = Module.cwrap("AdemAlgebra_basisElement_getBocksteins", 'number', ['pointer']);

// void AdemAlgebra_makeMonoAdmissible(AdemAlgebra *algebra, Vector *result, uint coeff, AdemBasisElement *monomial, int excess);
let cAdemAlgebra_makeMonoAdmissible = Module.cwrap("AdemAlgebra_makeMonoAdmissible", 'void', ['pointer', 'pointer', 'number', 'pointer', 'number']);



/**
 * milnor.h
 */

// Profile *Profile_construct(bool generic, uint q_part_length, uint * q_part, uint p_part_length, uint *p_part, bool truncated);
let cProfile_construct = cwrap("Profile_construct", 'pointer', ['bool', 'number', 'pointer', 'number', 'pointer', 'bool']);
let cProfile_free = cwrap("Profile_free", 'void', ['pointer']);

// Algebra *MilnorAlgebra_construct(uint p, bool generic, Profile *profile);
let cMilnorAlgebra_construct = Module.cwrap("MilnorAlgebra_construct", 'pointer', ['number', 'bool', 'pointer']);

// void MilnorAlgebra_free(MilnorAlgebra *);
let cMilnorAlgebra_free = Module.cwrap("MilnorAlgebra_free", 'void', ['pointer']);

// void MilnorAlgebra_generateBasis(Algebra *algebra, int max_degree);
let cMilnorAlgebra_generateBasis = Module.cwrap("MilnorAlgebra_generateBasis", 'void', ['pointer', 'number']);

// uint MilnorAlgebra_getDimension(Algebra *algebra, int degree, int excess);
let cMilnorAlgebra_getDimension = Module.cwrap("MilnorAlgebra_getDimension", 'number', ['pointer', 'number', 'number']);


/**
 * modules.h
 */
// uint Module_getDimension_function(Module *module, int degree);
let cModule_getDimension = Module.cwrap("Module_getDimension_function", 'number', ['pointer', 'number'])

// FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, char *name, int min_degree, int max_degree, uint *graded_dimension);
let cFiniteDimensionalModule_construct = Module.cwrap("FiniteDimensionalModule_construct", 'pointer', ['bool', 'pointer', 'pointer', 'number', 'number', 'pointer']);

// void FiniteDimensionalModule_free(FiniteDimensionalModule *module);
let cFiniteDimensionalModule_free = Module.cwrap("FiniteDimensionalModule_free", 'void', ['pointer']);

// void FiniteDimensionalModule_setAction(
//     FiniteDimensionalModule *module,
//     int operation_degree, uint operation_idx,
//     int input_degree, uint input_idx,
//     uint *output
// );
let cFiniteDimensionalModule_setAction = cwrap("FiniteDimensionalModule_setAction", 'void', 
    ['pointer', 'number' , 'number', 'number' , 'number', 'pointer']
);


// FreeModule *FreeModuleHomomorphism_getSource(FreeModuleHomomorphism *f);
let cFreeModuleHomomorphism_getSource = cwrap("FreeModuleHomomorphism_getSource", 'pointer', ['pointer']);

// Module *FreeModuleHomomorphism_getTarget(FreeModuleHomomorphism *f);
let cFreeModuleHomomorphism_getTarget = cwrap("FreeModuleHomomorphism_getTarget", 'pointer', ['pointer']);

// void FreeModuleHomomorphism_applyToGenerator(FreeModuleHomomorphism *f, Vector *result, uint coeff, int generator_degree, uint generator_index);
let cFreeModuleHomomorphism_applyToGenerator = cwrap("FreeModuleHomomorphism_applyToGenerator", 'void', ['pointer', 'pointer', 'number', 'number', 'number']);

// uint FreeModule_element_toJSONString(char *result, FreeModule *this, int degree, Vector *element)
let cFreeModule_element_toJSONString = cwrap("FreeModule_element_toJSONString", 'number', ['pointer', 'pointer', 'number', 'pointer']);

/**
 * resolution.h
 */

// Resolution *Resolution_construct(Module, max_degree, addClass, addStructline);
let cResolution_construct = Module.cwrap("Resolution_construct", 'pointer', ['pointer', 'number', 'pointer', 'pointer']);

// void Resolution_free(Resolution *res);
let cResolution_free = Module.cwrap("Resolution_free", 'void', ['pointer']);

// FreeModuleHomomorphism *Resolution_getDifferential(Resolution *resolution, uint homological_degree);
let cResolution_getDifferential = Module.cwrap("Resolution_getDifferential", 'pointer', ['pointer', 'number']);

// void Resolution_resolveThroughDegree(Resolution *res, uint degree);
let cResolution_resolveThroughDegree = Module.cwrap("Resolution_resolveThroughDegree", 'void', ['pointer', 'number']);