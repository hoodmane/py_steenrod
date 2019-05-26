cmake_minimum_required(VERSION 3.9)
project(CSteenrod VERSION 1.0.1 DESCRIPTION "C Steenrod Modules")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(DEBUG_INFO_COMPILE_FLAGS "-g")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${DEBUG_INFO_COMPILE_FLAGS}")
SET(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -lprofiler")


include(GNUInstallDirs)
add_executable(csteen
        include/algebra.h
        include/combinatorics.h
        include/milnor.h
        src/khash.h
        src/combinatorics.c
        src/combinatorics.h
        src/FpVector.c
        src/FpVector.h
        src/algebra.h
        src/modules.h
        src/modules.c
        src/resolution.h
        src/resolution.c
        src/milnor.c
        src/milnor.h
        src/milnor_datatypes.c
        src/milnor_private.h
        src/modules.c
    )
add_library(CSteenrod SHARED
        include/algebra.h
        include/combinatorics.h
        include/milnor.h
        src/khash.h
        src/combinatorics.c
        src/combinatorics.h
        src/FpVector.c
        src/FpVector.h
        src/algebra.h
        src/modules.h
        src/modules.c
        src/resolution.h
        src/resolution.c
        src/milnor.c
        src/milnor.h
        src/milnor_datatypes.c
        src/milnor_private.h
        src/modules.c
    )

set_target_properties(CSteenrod PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)
target_include_directories(CSteenrod PUBLIC include PRIVATE src)

configure_file(CSteenrod.pc.in CSteenrod.pc @ONLY)
target_include_directories(CSteenrod PRIVATE .)
install(TARGETS CSteenrod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
