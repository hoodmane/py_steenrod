// void initializePrime(uint p);
let cinitializePrime = cwrap("initializePrime", 'void', ['number']);

// Vector *Vector_construct(uint p, uint dimension, uint offset);
let cVector_construct = cwrap("Vector_construct", 'pointer', ['number', 'number', 'number']);

// void Vector_pack(Vector *target, uint *source);
let cVector_pack = cwrap("Vector_pack", 'void', ['pointer', 'pointer']);

// void Vector_unpack(uint * target, Vector * source);
let cVector_unpack = cwrap("Vector_unpack", 'void', ['pointer', 'pointer']);

// Profile *Profile_construct(bool generic, uint q_part_length, uint * q_part, uint p_part_length, uint *p_part, bool truncated);
let cProfile_construct = cwrap("Profile_construct", 'pointer', ['bool', 'number', 'pointer', 'number', 'pointer', 'bool']);
let cProfile_free = cwrap("Profile_free", 'void', ['pointer']);

// Algebra *MilnorAlgebra_construct(uint p, bool generic, Profile *profile);
let cMilnorAlgebra_construct = Module.cwrap("MilnorAlgebra_construct", 'pointer', ['number', 'bool', 'pointer']);

// void MilnorAlgebra_free(MilnorAlgebra *);
let cMilnorAlgebra_free = Module.cwrap("MilnorAlgebra_free", 'void', ['pointer']);

// bool MilnorAlgebra_generateBasis(Algebra *algebra, uint max_degree);
let cMilnorAlgebra_generateBasis = Module.cwrap("MilnorAlgebra_generateBasis", 'bool', ['pointer', 'number']);

// uint MilnorAlgebra_getDimension(Algebra *algebra, uint degree);
let cMilnorAlgebra_getDimension = Module.cwrap("MilnorAlgebra_getDimension", 'number', ['pointer', 'number']);

// FiniteDimensionalModule *FiniteDimensionalModule_construct(Algebra *algebra, uint max_generator_degree, uint *graded_dimension);
let cFiniteDimensionalModule_construct = Module.cwrap("FiniteDimensionalModule_construct", 'pointer', ['bool', 'pointer', 'number', 'pointer']);

// void FiniteDimensionalModule_free(FiniteDimensionalModule *module);
let cFiniteDimensionalModule_free = Module.cwrap("FiniteDimensionalModule_free", 'void', ['pointer']);

// void FiniteDimensionalModule_setAction(
//     FiniteDimensionalModule *module,
//     uint operation_degree, uint operation_idx,
//     uint input_degree, uint input_idx,
//     uint *output
// );
let cFiniteDimensionalModule_setAction = cwrap("FiniteDimensionalModule_setAction", 'void', 
    ['pointer', 'number' , 'number', 'number' , 'number', 'pointer']
);

let cdoResolution = Module.cwrap("doResolution", 'void', ['number', 'pointer', 'pointer']);

//Resolution *Resolution_construct(Module, max_degree, addClass, addStructline)
let cResolution_construct = Module.cwrap("Resolution_construct", 'pointer', ['pointer', 'number', 'pointer', 'pointer']);

// void Resolution_resolveThroughDegree(Resolution *res, uint degree);
let cResolution_resolveThroughDegree = Module.cwrap("Resolution_resolveThroughDegree", 'void', ['pointer', 'number']);