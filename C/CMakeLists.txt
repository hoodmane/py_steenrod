cmake_minimum_required(VERSION 3.9)
project(CSteenrod VERSION 1.0.1 DESCRIPTION "C Steenrod Modules")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(DEBUG_INFO_COMPILE_FLAGS "-g")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${DEBUG_INFO_COMPILE_FLAGS}")


include(GNUInstallDirs)
add_executable(csteen src/milnor.c src/milnor_datatypes.c src/combinatorics.c src/algebra.c)
add_library(CSteenrod SHARED src/milnor.c src/milnor_datatypes.c src/combinatorics.c src/algebra.c)
set_target_properties(CSteenrod PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER src/combinatorics.h src/milnor.h src/milnor_datatypes.h)
configure_file(CSteenrod.pc.in CSteenrod.pc @ONLY)
target_include_directories(CSteenrod PRIVATE .)
install(TARGETS CSteenrod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
