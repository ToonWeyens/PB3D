cmake_minimum_required(VERSION 3.18)

project(PB3D
    VERSION 2.48
    LANGUAGES Fortran C CXX
    DESCRIPTION "Peeling Ballooning in 3D"
)

# Require out-of-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Please create a build directory.")
endif()

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Options
option(PB3D_ENABLE_DEBUG "Enable debug mode (ldebug preprocessor flag)" OFF)
option(PB3D_ENABLE_INFINIBAND "Enable InfiniBand support (lIB preprocessor flag)" OFF)
option(PB3D_BUILD_EXECUTABLES "Build the PB3D and POST executables (requires the full dependency stack: PETSc, SLEPc, HDF5, NetCDF, PSPLINE, LIBSTELL, ...)" ON)

# CTest integration; provides the BUILD_TESTING option (default ON).
# With PB3D_BUILD_EXECUTABLES=OFF, only the dependency-light core library and
# its unit tests are built, which requires no external libraries at all.
include(CTest)

# Build type defaults to Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#------------------------------------------------------------------------------
# Detect compiler
#------------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(PB3D_COMPILER_GNU TRUE)
    set(PB3D_COMPILER_INTEL FALSE)
    message(STATUS "Detected GNU Fortran compiler")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(PB3D_COMPILER_GNU FALSE)
    set(PB3D_COMPILER_INTEL TRUE)
    message(STATUS "Detected Intel Fortran compiler")
else()
    message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Assuming GNU-compatible.")
    set(PB3D_COMPILER_GNU TRUE)
    set(PB3D_COMPILER_INTEL FALSE)
endif()

#------------------------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------------------------
# Enable preprocessing for all Fortran files
set(CMAKE_Fortran_PREPROCESS ON)

# GNU Fortran flags
# Note: -fimplicit-none is intentionally not global to support legacy libraries like dfftpack
if(PB3D_COMPILER_GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fbacktrace -fno-omit-frame-pointer -fallow-argument-mismatch")
    # no -pedantic: it re-promotes the argument mismatches in vendored F77
    # (Libraries/dfft.f) to hard errors despite -fallow-argument-mismatch
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -Wextra -finit-real=snan -ffpe-trap=invalid -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
endif()

# Intel Fortran flags
if(PB3D_COMPILER_INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -implicitnone -traceback -fno-omit-frame-pointer -heap-arrays 100 -recursive")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -ftrapuv -check bounds -check uninit -diag-disable 6536 -diag-disable 6843")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -diag-disable 6536 -diag-disable 6843")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -diag-disable 6536 -diag-disable 6843")
endif()

# Preprocessor definitions
if(PB3D_ENABLE_DEBUG)
    add_compile_definitions(ldebug)
    message(STATUS "Debug mode enabled (ldebug)")
endif()

if(PB3D_ENABLE_INFINIBAND)
    add_compile_definitions(lIB)
    message(STATUS "InfiniBand support enabled (lIB)")
endif()

# Compiler-specific preprocessor defines
if(PB3D_COMPILER_INTEL)
    add_compile_definitions(lwith_intel)
else()
    add_compile_definitions(lwith_gnu)
endif()

#------------------------------------------------------------------------------
# Find dependencies
#------------------------------------------------------------------------------

# MPI: needed by the full build always, and by the core library only in debug
# mode (messages.f90 uses MPI for memory monitoring under ldebug)
if(PB3D_BUILD_EXECUTABLES OR PB3D_ENABLE_DEBUG)
    find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
    message(STATUS "Found MPI: ${MPI_Fortran_COMPILER}")
endif()

if(PB3D_BUILD_EXECUTABLES)

# PETSc (uses its own CMake config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
message(STATUS "Found PETSc: ${PETSC_LINK_LIBRARIES}")

# SLEPc (uses its own CMake config)
pkg_check_modules(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
message(STATUS "Found SLEPc: ${SLEPC_LINK_LIBRARIES}")

# HDF5
find_package(HDF5 REQUIRED COMPONENTS Fortran)
message(STATUS "Found HDF5: ${HDF5_Fortran_LIBRARIES}")

# NetCDF-Fortran
pkg_check_modules(NETCDF_FORTRAN REQUIRED IMPORTED_TARGET netcdf-fortran)
message(STATUS "Found NetCDF-Fortran: ${NETCDF_FORTRAN_LINK_LIBRARIES}")

# ScaLAPACK
find_package(PkgConfig REQUIRED)
pkg_check_modules(SCALAPACK IMPORTED_TARGET scalapack)
if(NOT SCALAPACK_FOUND)
    # Try to find ScaLAPACK library directly
    find_library(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi)
    if(SCALAPACK_LIBRARY)
        set(SCALAPACK_FOUND TRUE)
        set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
        message(STATUS "Found ScaLAPACK: ${SCALAPACK_LIBRARY}")
    else()
        message(WARNING "ScaLAPACK not found. Some features may not work.")
    endif()
endif()

# STRUMPACK-Dense 1.1.1 (required for vacuum module)
# This is the old dense-only version with Fortran interface
# Set STRUMPACK_DIR to the STRUMPACK-Dense-1.1.1 directory
if(NOT STRUMPACK_DIR)
    set(STRUMPACK_DIR "$ENV{HOME}/Code/STRUMPACK-Dense-1.1.1")
endif()
# Use NO_DEFAULT_PATH to avoid finding the newer STRUMPACK from Spack
find_library(STRUMPACK_LIBRARY NAMES strumpack
    PATHS ${STRUMPACK_DIR}/lib ${STRUMPACK_DIR}/src
    NO_DEFAULT_PATH)
find_path(STRUMPACK_INCLUDE_DIR NAMES strumpackdensepackage.mod
    PATHS ${STRUMPACK_DIR}/inc ${STRUMPACK_DIR}/examples
    NO_DEFAULT_PATH)
if(STRUMPACK_LIBRARY AND STRUMPACK_INCLUDE_DIR)
    set(STRUMPACK_FOUND TRUE)
    set(STRUMPACK_LIBRARIES ${STRUMPACK_LIBRARY})
    set(STRUMPACK_INCLUDE_DIRS ${STRUMPACK_INCLUDE_DIR})
    message(STATUS "Found STRUMPACK-Dense: ${STRUMPACK_LIBRARY}")
    message(STATUS "  Include dir: ${STRUMPACK_INCLUDE_DIR}")
else()
    set(STRUMPACK_FOUND FALSE)
    message(STATUS "STRUMPACK-Dense not found - the vacuum module will solve its \
boundary-element system with ScaLAPACK (pdgesv) instead of STRUMPACK's \
compressed solver. Set STRUMPACK_DIR to a STRUMPACK-Dense-1.1.1 build to \
use the latter.")
endif()

# PSPLINE (custom find module)
find_package(PSPLINE REQUIRED)
message(STATUS "Found PSPLINE: ${PSPLINE_LIBRARIES}")

# LIBSTELL (custom find module)
find_package(LIBSTELL REQUIRED)
message(STATUS "Found LIBSTELL: ${LIBSTELL_LIBRARIES}")

# LAPACK and BLAS (usually pulled in by other dependencies, but ensure they're available)
find_package(LAPACK)
find_package(BLAS)

endif() # PB3D_BUILD_EXECUTABLES

#------------------------------------------------------------------------------
# Internal libraries
#------------------------------------------------------------------------------

# libdfftpack - from Libraries/dfft.f
add_library(dfftpack STATIC
    ${CMAKE_SOURCE_DIR}/Libraries/dfft.f
)
target_compile_options(dfftpack PRIVATE -O2 -funroll-loops -fexpensive-optimizations -fallow-argument-mismatch)

# libfoul - from Libraries/foul.f90
add_library(foul STATIC
    ${CMAKE_SOURCE_DIR}/Libraries/foul.f90
)
target_compile_options(foul PRIVATE -O2 -w)
set_target_properties(foul PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/foul)
target_include_directories(foul PUBLIC ${CMAKE_BINARY_DIR}/mod/foul)

#------------------------------------------------------------------------------
# Dependency-light core library
#
# The leaf modules of the dependency tree: everything here compiles with just
# a Fortran compiler (plus MPI in debug mode), without PETSc/SLEPc/HDF5/
# NetCDF/PSPLINE/LIBSTELL. The unit tests in tests/ build against this
# library, so they remain available on machines without the full stack.
#------------------------------------------------------------------------------
set(PB3D_CORE_SOURCES
    Modules/num_vars.f90
    Modules/str_utilities.f90
    Modules/messages.f90
    Modules/files_utilities.f90
    Modules/dtorh.f90
)

# num_utilities (pure numerics: integration, finite-difference weights, index
# helpers, zero finders, ...) needs LAPACK for its 0-D determinant/inverse and
# polynomial extrapolation; since the spline wrappers were split off into
# spline_utilities it has no other external dependency. num_ops additionally
# uses output_ops (HDF5) under ldebug, so it only joins the core in non-debug
# configurations.
find_package(LAPACK QUIET)
set(PB3D_CORE_HAS_NUM FALSE)
if(LAPACK_FOUND)
    list(APPEND PB3D_CORE_SOURCES Modules/num_utilities.f90)
    set(PB3D_CORE_HAS_NUM TRUE)
    if(NOT PB3D_ENABLE_DEBUG)
        list(APPEND PB3D_CORE_SOURCES Modules/num_ops.f90)
        set(PB3D_CORE_HAS_NUM_OPS TRUE)
    endif()
else()
    message(STATUS "LAPACK not found - num_utilities/num_ops unit tests disabled")
endif()

add_library(pb3d_core OBJECT ${PB3D_CORE_SOURCES})
set_target_properties(pb3d_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/core)
target_include_directories(pb3d_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/mod/core
)
target_link_libraries(pb3d_core PUBLIC foul)
if(LAPACK_FOUND)
    target_link_libraries(pb3d_core PUBLIC ${LAPACK_LIBRARIES})
endif()
if(PB3D_ENABLE_DEBUG)
    # messages.f90 uses MPI under ldebug (memory monitoring)
    target_link_libraries(pb3d_core PUBLIC MPI::MPI_Fortran)
endif()

if(PB3D_BUILD_EXECUTABLES)

#------------------------------------------------------------------------------
# Module source files (order matters for Fortran module dependencies)
#------------------------------------------------------------------------------
set(PB3D_MODULE_SOURCES
    # Base modules (no dependencies on other PB3D modules)
    Modules/num_vars.f90
    Modules/str_utilities.f90

    # Messages and utilities
    Modules/messages.f90
    Modules/num_utilities.f90
    Modules/spline_utilities.f90
    Modules/num_ops.f90
    Modules/files_utilities.f90

    # MPI modules
    Modules/MPI_vars.f90
    Modules/MPI_utilities.f90
    Modules/MPI_ops.f90

    # HDF5 modules
    Modules/HDF5_vars.f90
    Modules/HDF5_utilities.f90
    Modules/HDF5_ops.f90

    # Grid modules
    Modules/grid_vars.f90
    Modules/eq_vars.f90
    Modules/grid_utilities.f90
    Modules/grid_ops.f90

    # Equilibrium modules
    Modules/eq_utilities.f90
    Modules/eq_ops.f90

    # VMEC/HELENA modules
    Modules/VMEC_vars.f90
    Modules/VMEC_utilities.f90
    Modules/VMEC_ops.f90
    Modules/HELENA_vars.f90
    Modules/HELENA_ops.f90

    # Perturbation modules
    Modules/X_vars.f90
    Modules/X_utilities.f90
    Modules/X_ops.f90

    # Solution modules
    Modules/sol_vars.f90
    Modules/sol_utilities.f90
    Modules/sol_ops.f90

    # Vacuum modules
    Modules/vac_vars.f90
    Modules/vac_utilities.f90
    Modules/vac_ops.f90
    Modules/dtorh.f90

    # Richardson extrapolation
    Modules/rich_vars.f90
    Modules/rich_ops.f90

    # SLEPc interface
    Modules/SLEPC_utilities.f90
    Modules/SLEPC_ops.f90

    # I/O modules
    Modules/input_utilities.f90
    Modules/input_ops.f90
    Modules/output_ops.f90
    Modules/files_ops.f90

    # PB3D modules
    Modules/PB3D_utilities.f90
    Modules/PB3D_ops.f90

    # Driver modules
    Modules/driver_eq.f90
    Modules/driver_X.f90
    Modules/driver_sol.f90
    Modules/driver_POST.f90

    # Test module
    Modules/test.f90
)

#------------------------------------------------------------------------------
# PB3D object library (shared between PB3D and POST)
#------------------------------------------------------------------------------
add_library(pb3d_modules OBJECT ${PB3D_MODULE_SOURCES})

# Keep the .mod files of the full module library separate from those of
# pb3d_core (the core sources are compiled into both targets)
set_target_properties(pb3d_modules PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/full)

# Include directories
target_include_directories(pb3d_modules PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/mod/full
    ${PSPLINE_INCLUDE_DIRS}
    ${LIBSTELL_INCLUDE_DIRS}
    ${HDF5_Fortran_INCLUDE_DIRS}
)

# Link dependencies (for module interface)
target_link_libraries(pb3d_modules PUBLIC
    foul
    MPI::MPI_Fortran
    PkgConfig::PETSC
    PkgConfig::SLEPC
    ${HDF5_Fortran_LIBRARIES}
    PkgConfig::NETCDF_FORTRAN
)

if(STRUMPACK_FOUND)
    target_include_directories(pb3d_modules PUBLIC ${STRUMPACK_INCLUDE_DIRS})
    target_compile_definitions(pb3d_modules PRIVATE PB3D_WITH_STRUMPACK)
endif()

#------------------------------------------------------------------------------
# PB3D executable
#------------------------------------------------------------------------------
add_executable(PB3D
    PB3D.f90
    $<TARGET_OBJECTS:pb3d_modules>
)

target_include_directories(PB3D PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/mod/full
    ${PSPLINE_INCLUDE_DIRS}
    ${LIBSTELL_INCLUDE_DIRS}
    ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(PB3D PRIVATE
    dfftpack
    foul
    ${PSPLINE_LIBRARIES}
    ${LIBSTELL_LIBRARIES}
    MPI::MPI_Fortran
    MPI::MPI_CXX
    PkgConfig::PETSC
    PkgConfig::SLEPC
    ${HDF5_Fortran_LIBRARIES}
    PkgConfig::NETCDF_FORTRAN
)

if(SCALAPACK_FOUND)
    if(TARGET PkgConfig::SCALAPACK)
        target_link_libraries(PB3D PRIVATE PkgConfig::SCALAPACK)
    else()
        target_link_libraries(PB3D PRIVATE ${SCALAPACK_LIBRARIES})
    endif()
endif()

if(STRUMPACK_FOUND)
    target_include_directories(PB3D PRIVATE ${STRUMPACK_INCLUDE_DIRS})
    target_link_libraries(PB3D PRIVATE ${STRUMPACK_LIBRARIES})
endif()

if(LAPACK_FOUND)
    target_link_libraries(PB3D PRIVATE ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
    target_link_libraries(PB3D PRIVATE ${BLAS_LIBRARIES})
endif()

# Link standard C++ and math libraries
target_link_libraries(PB3D PRIVATE m stdc++)

#------------------------------------------------------------------------------
# POST executable
#------------------------------------------------------------------------------
add_executable(POST
    POST.f90
    $<TARGET_OBJECTS:pb3d_modules>
)

target_include_directories(POST PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/mod/full
    ${PSPLINE_INCLUDE_DIRS}
    ${LIBSTELL_INCLUDE_DIRS}
    ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(POST PRIVATE
    dfftpack
    foul
    ${PSPLINE_LIBRARIES}
    ${LIBSTELL_LIBRARIES}
    MPI::MPI_Fortran
    MPI::MPI_CXX
    PkgConfig::PETSC
    PkgConfig::SLEPC
    ${HDF5_Fortran_LIBRARIES}
    PkgConfig::NETCDF_FORTRAN
)

if(SCALAPACK_FOUND)
    if(TARGET PkgConfig::SCALAPACK)
        target_link_libraries(POST PRIVATE PkgConfig::SCALAPACK)
    else()
        target_link_libraries(POST PRIVATE ${SCALAPACK_LIBRARIES})
    endif()
endif()

if(STRUMPACK_FOUND)
    target_include_directories(POST PRIVATE ${STRUMPACK_INCLUDE_DIRS})
    target_link_libraries(POST PRIVATE ${STRUMPACK_LIBRARIES})
endif()

if(LAPACK_FOUND)
    target_link_libraries(POST PRIVATE ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
    target_link_libraries(POST PRIVATE ${BLAS_LIBRARIES})
endif()

# Link standard C++ and math libraries
target_link_libraries(POST PRIVATE m stdc++)

#------------------------------------------------------------------------------
# Installation
#------------------------------------------------------------------------------
install(TARGETS PB3D POST
    RUNTIME DESTINATION bin
)

endif() # PB3D_BUILD_EXECUTABLES

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------
message(STATUS "")
message(STATUS "PB3D Configuration Summary")
message(STATUS "==========================")
message(STATUS "  Version:        ${PROJECT_VERSION}")
message(STATUS "  Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran:        ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Compiler ID:    ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  Debug mode:     ${PB3D_ENABLE_DEBUG}")
message(STATUS "  InfiniBand:     ${PB3D_ENABLE_INFINIBAND}")
message(STATUS "  Executables:    ${PB3D_BUILD_EXECUTABLES}")
message(STATUS "  Testing:        ${BUILD_TESTING}")
message(STATUS "")
if(PB3D_BUILD_EXECUTABLES)
    message(STATUS "Dependencies:")
    message(STATUS "  MPI:            ${MPI_Fortran_COMPILER}")
    message(STATUS "  PETSc:          ${PETSC_VERSION}")
    message(STATUS "  SLEPc:          ${SLEPC_VERSION}")
    message(STATUS "  HDF5:           ${HDF5_VERSION}")
    message(STATUS "  NetCDF-Fortran: ${NETCDF_FORTRAN_VERSION}")
    message(STATUS "  ScaLAPACK:      ${SCALAPACK_FOUND}")
    message(STATUS "  STRUMPACK:      ${STRUMPACK_FOUND}")
    message(STATUS "  PSPLINE:        ${PSPLINE_FOUND}")
    message(STATUS "  LIBSTELL:       ${LIBSTELL_FOUND}")
    message(STATUS "")
endif()
