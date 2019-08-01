cmake_minimum_required(VERSION 3.9)
project(CSteenrod VERSION 1.0.1 DESCRIPTION "C Steenrod Modules")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(DEBUG_COMPILE_FLAGS -Wno-unused-variable -Wall -Wextra -Wpedantic -Werror -O0 -g3)
SET(OPTIMIZE_COMPILE_FLAGS -O3 -DNDEBUG) #-pg)
SET(CMAKE_CXX_FLAGS "")
#SET(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -pg") #-lprofiler
#SET(TARGET_COMPILE_OPTIONS  "${DEBUG_COMPILE_FLAGS}")
SET(TARGET_COMPILE_OPTIONS "${OPTIMIZE_COMPILE_FLAGS}")


include(GNUInstallDirs)
add_executable(csteen
        src/khash.h
        src/parson.h
        src/parson.c

        src/combinatorics.h
        src/combinatorics.c
        src/FpVector.h
        src/FpVector.c
        src/Matrix.h
        src/Matrix.c

        src/Algebra.h
        src/Algebra.c
        src/MilnorAlgebra.h
        src/milnor_private.h
        src/MilnorAlgebra.c
        src/milnor_io.c
        src/AdemAlgebra.h
        src/AdemAlgebra.c
        src/adem_io.c

        src/Module.h
        src/Module.c  
        src/FiniteDimensionalModule.h
        src/FiniteDimensionalModule.c
        src/FreeModule.h
        src/FreeModule.c
        src/FreeModuleHomomorphism.h
        src/FreeModuleHomomorphism.c
        src/Resolution.h
        src/Resolution.c
        src/ResolutionHomomorphism.h
        src/ResolutionHomomorphism.c     
        src/main.h
        src/main.c        
    )
add_library(CSteenrod SHARED
        src/khash.h
        src/parson.h
        src/parson.c

        src/combinatorics.h
        src/combinatorics.c
        src/FpVector.h
        src/FpVector.c
        src/Matrix.h
        src/Matrix.c

        src/Algebra.h
        src/Algebra.c
        src/MilnorAlgebra.h
        src/milnor_private.h
        src/MilnorAlgebra.c
        src/milnor_io.c
        src/AdemAlgebra.h
        src/AdemAlgebra.c
        src/adem_io.c

        src/Module.h
        src/Module.c  
        src/FiniteDimensionalModule.h
        src/FiniteDimensionalModule.c
        src/FreeModule.h
        src/FreeModule.c
        src/FreeModuleHomomorphism.h
        src/FreeModuleHomomorphism.c
        src/Resolution.h
        src/Resolution.c
        src/ResolutionHomomorphism.h
        src/ResolutionHomomorphism.c
        src/main.h
        src/main.c
    )

set_target_properties(CSteenrod PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)
target_include_directories(CSteenrod PUBLIC include PRIVATE src)
target_compile_options(CSteenrod PUBLIC ${TARGET_COMPILE_OPTIONS})
target_compile_options(csteen PUBLIC ${TARGET_COMPILE_OPTIONS})
target_compile_options(csteen PUBLIC -O0 -g3)
#target_link_options(csteen PUBLIC "-pg")

configure_file(CSteenrod.pc.in CSteenrod.pc @ONLY)
target_include_directories(CSteenrod PRIVATE .)
install(TARGETS CSteenrod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
