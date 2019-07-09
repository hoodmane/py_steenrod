cmake_minimum_required(VERSION 3.9)
project(CSteenrod VERSION 1.0.1 DESCRIPTION "C Steenrod Modules")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(DEBUG_INFO_COMPILE_FLAGS "-g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_INFO_COMPILE_FLAGS}")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -s TOTAL_MEMORY=50331648") # 3 * 2^23

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g4 -O0 --source-map-base http://localhost:8080/Cemsdk")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -O3")

SET(EXTRA_EXPORTED_RUNTIME_METHODS 
    "['ccall', 'cwrap', 'addFunction', 'print', 'getValue', 'setValue']"
)
SET(EXPORTED_FUNCTIONS "'_initializePrime'")
SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_Vector_construct', '_Vector_pack', '_Vector_unpack', '_Vector_print', '_Vector_free'")

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_algebra_computeBasis_function', '_algebra_getDimension_function', \
    '_algebra_multiplyBasisElements_function', '_algebra_basisElementToString_function'")

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_MilnorAlgebra_construct', '_MilnorAlgebra_free', \
    '_MilnorAlgebra_generateBasis', '_MilnorAlgebra_getDimension', \
    '_Profile_construct', '_Profile_free'"
)

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_AdemAlgebra_construct', '_AdemAlgebra_free', \
    '_AdemAlgebra_generateBasis', '_AdemAlgebra_getDimension', \
    '_AdemAlgebra_basisElement_fromIndex', '_AdemAlgebra_basisElement_toIndex', \
    '_AdemAlgebra_multiply', '_AdemAlgebra_makeMonoAdmissible', \
    '_AdemAlgebra_element_toString', \
    '_AdemAlgebra_basisElement_construct', '_AdemAlgebra_basisElement_free', \
    '_AdemAlgebra_basisElement_getPlength', '_AdemAlgebra_basisElement_getPs', \
    '_AdemAlgebra_basisElement_getBocksteins', '_AdemAlgebra_basisElement_fromIndex'"
)

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_Module_getDimension_function'"
)

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_FiniteDimensionalModule_construct', '_FiniteDimensionalModule_free', \
    '_FiniteDimensionalModule_setAction'"
)

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_FreeModule_getDimension', '_FreeModuleHomomorphism_applyToGenerator', \
    '_FreeModuleHomomorphism_getSource', '_FreeModuleHomomorphism_getTarget', \
    '_FreeModule_element_toJSONString'"
)

SET(EXPORTED_FUNCTIONS "${EXPORTED_FUNCTIONS},\
    '_Resolution_construct', '_Resolution_resolveThroughDegree', \
    '_Resolution_getDifferential'"
)

SET(EMCC_EXPORT_FLAGS "-s EXTRA_EXPORTED_RUNTIME_METHODS=\"${EXTRA_EXPORTED_RUNTIME_METHODS}\"")
SET(EMCC_EXPORT_FLAGS "${EMCC_EXPORT_FLAGS} -s EXPORTED_FUNCTIONS=\"[${EXPORTED_FUNCTIONS}]\"")
SET(EMCC_EXPORT_FLAGS "${EMCC_EXPORT_FLAGS} -s RESERVED_FUNCTION_POINTERS=20")
SET(EMCC_EXPORT_FLAGS "${EMCC_EXPORT_FLAGS} -s TOTAL_MEMORY=167772160")


include(GNUInstallDirs)
add_executable(CSteenrod
        ../C/src/khash.h
        ../C/src/combinatorics.h
        ../C/src/combinatorics.c
        ../C/src/FpVector.h
        ../C/src/FpVector.c
        ../C/src/algebra.h
        ../C/src/algebra.c
        ../C/src/Module.h
        ../C/src/Module.c  
        ../C/src/FiniteDimensionalModule.h
        ../C/src/FiniteDimensionalModule.c
        ../C/src/FreeModule.h
        ../C/src/FreeModule.c
        ../C/src/FreeModuleHomomorphism.h
        ../C/src/FreeModuleHomomorphism.c
        ../C/src/resolution.h
        ../C/src/resolution.c
        ../C/src/milnor.h
        ../C/src/milnor_private.h
        ../C/src/milnor.c
        ../C/src/milnor_io.c
        ../C/src/adem.h
        ../C/src/adem.c
        ../C/src/adem_io.c
    )

if (${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
    set_target_properties(CSteenrod PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -s ASSERTIONS=1 -g4 --source-map-base http://localhost:8080/Cemsdk")
    set_target_properties(CSteenrod PROPERTIES LINK_FLAGS "${LINK_FLAGS} -g4 ${EMCC_EXPORT_FLAGS}")
endif ()

configure_file(../C/CSteenrod.pc.in ../C/CSteenrod.pc @ONLY)
