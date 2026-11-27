cmake_minimum_required(VERSION 3.18)

project(PB3D
    VERSION 2.57
    LANGUAGES Fortran C
    DESCRIPTION "Peeling Ballooning in 3D"
)

# Require out-of-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Please create a build directory.")
endif()

# Add cmake modules path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# The Makefile generator finds the dependencies of a Fortran source with
# CMake's own preprocessor, which follows #include into external headers too.
# CMake 3.31 never finishes on PETSc 3.25's Fortran headers, where
# petscsys.h and petscmath.h include each other under #if guards: on
# Snellius, `cmake --build` scanned for a quarter of an hour at 13 GB without
# compiling anything. With the scan limited to the files of this project,
# CMake 3.31 builds normally. The price: a change in an external header does
# not by itself rebuild PB3D, so clean the build directory after reinstalling
# PETSc in place. The Ninja generator ignores this variable.
set(CMAKE_DEPENDS_IN_PROJECT_ONLY ON)

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
    # gfortran 14.1 to 14.3 and 15.1 to 15.2 miscompile the bounds check of a
    # subscript that calls a function with an array-constructor argument
    # inside an actual argument, as in eq_2%T_EF(:,:,:,c([1,1],.false.),...)
    # throughout eq_ops: the check calls c a second time without setting up
    # the descriptor of the constructor, so it reads a garbage index or
    # segfaults inside c (GCC PR fortran/119273, fixed in 14.4, 15.3 and 16.1;
    # standalone reproducer: tests/reference/gfortran14_bounds_repro.f90).
    # Those versions get a Debug build without -fcheck=bounds instead of none
    # at all.
    if((CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 14 AND
        CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 14.4) OR
       (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 15 AND
        CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 15.3))
        string(REPLACE "-fcheck=all" "-fcheck=all,no-bounds"
               CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
        message(STATUS "gfortran ${CMAKE_Fortran_COMPILER_VERSION}: Debug build without -fcheck=bounds (compiler bug, see CMakeLists.txt)")
    endif()
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

# macOS: PB3D warns when the allocator keeps freed memory (input_ops.f90)
if(APPLE)
    add_compile_definitions(lmacos)
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
    find_package(MPI REQUIRED COMPONENTS Fortran C)
    message(STATUS "Found MPI: ${MPI_Fortran_COMPILER}")
endif()

if(PB3D_BUILD_EXECUTABLES)

# PETSc and SLEPc, through pkg-config: PETSc.pc and SLEPc.pc in
# <prefix>/lib/pkgconfig must be on PKG_CONFIG_PATH (a Spack environment and
# scripts/snellius/env.sh put them there). cmake/PB3DCheckPETSc.cmake then
# checks complex scalars, the versions and the distributed LU.
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_check_modules(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
message(STATUS "Found PETSc ${PETSC_VERSION}: ${PETSC_PREFIX}")
message(STATUS "Found SLEPc ${SLEPC_VERSION}: ${SLEPC_PREFIX}")
include(PB3DCheckPETSc)

# HDF5 with the Fortran interface, built for MPI: PB3D writes its output
# collectively from all processes (H5Pset_fapl_mpio_f in HDF5_ops.f90)
find_package(HDF5 REQUIRED COMPONENTS Fortran)
message(STATUS "Found HDF5 ${HDF5_VERSION}: ${HDF5_Fortran_LIBRARIES}")
if(DEFINED HDF5_IS_PARALLEL AND NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "The HDF5 found is serial. PB3D needs HDF5 built \
with MPI (Spack: hdf5 +mpi +fortran, Ubuntu: libhdf5-openmpi-dev with \
-DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/openmpi).")
endif()

# NetCDF-Fortran: LIBSTELL reads the wout files of VMEC with it
pkg_check_modules(NETCDF_FORTRAN REQUIRED IMPORTED_TARGET netcdf-fortran)
message(STATUS "Found NetCDF-Fortran ${NETCDF_FORTRAN_VERSION}: ${NETCDF_FORTRAN_PREFIX}")

# ScaLAPACK with BLACS, required (cmake/FindSCALAPACK.cmake)
find_package(SCALAPACK REQUIRED)

# PSPLINE and LIBSTELL, built outside any package manager
# (cmake/FindPSPLINE.cmake, cmake/FindLIBSTELL.cmake)
find_package(PSPLINE REQUIRED)
find_package(LIBSTELL REQUIRED)

# LAPACK and BLAS (usually pulled in by other dependencies, but ensure they're available)
find_package(LAPACK)
find_package(BLAS)

#------------------------------------------------------------------------------
# All external libraries of the full build in one interface target, so that
# PB3D, POST and the full-stack tests (tests/CMakeLists.txt) link the same
# list in the same order
#------------------------------------------------------------------------------
add_library(pb3d_deps INTERFACE)
target_include_directories(pb3d_deps INTERFACE
    ${CMAKE_SOURCE_DIR}/include
    ${PSPLINE_INCLUDE_DIRS}
    ${LIBSTELL_INCLUDE_DIRS}
    ${HDF5_Fortran_INCLUDE_DIRS}
)
target_link_libraries(pb3d_deps INTERFACE
    ${PSPLINE_LIBRARIES}
    ${LIBSTELL_LIBRARIES}
    MPI::MPI_Fortran
    PkgConfig::PETSC
    PkgConfig::SLEPC
    ${HDF5_Fortran_LIBRARIES}
    PkgConfig::NETCDF_FORTRAN
    SCALAPACK::SCALAPACK
)
if(LAPACK_FOUND)
    target_link_libraries(pb3d_deps INTERFACE ${LAPACK_LIBRARIES})
endif()
if(BLAS_FOUND)
    target_link_libraries(pb3d_deps INTERFACE ${BLAS_LIBRARIES})
endif()
# libm
target_link_libraries(pb3d_deps INTERFACE m)

endif() # PB3D_BUILD_EXECUTABLES

#------------------------------------------------------------------------------
# Internal libraries
#------------------------------------------------------------------------------

# libdfftpack - from Libraries/dfft.f
add_library(dfftpack STATIC
    ${CMAKE_SOURCE_DIR}/Libraries/dfft.f
)
# -w: vendored code is not ours to lint (same policy as foul below and
# fortitude.toml); without it, dfft.f alone contributes ~157 warnings under
# -Wall -Wextra that the CI warning counter cannot reliably attribute,
# because gfortran prints the file location on a separate line
# -fcheck=no-bounds: dfft.f declares its work arrays with the F77 idiom
# DIMENSION WSAVE(1), so the Debug flag -fcheck=all stops every FFT call with
# a false bounds error (seen in the free-boundary regression deck, issue #4)
target_compile_options(dfftpack PRIVATE -O2 -funroll-loops -fexpensive-optimizations -fallow-argument-mismatch -w -fcheck=no-bounds)

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
    # light type modules and the pure index logic built on them
    Modules/flat_vars.f90
    Modules/grid_vars.f90
    Modules/X_vars.f90
    Modules/PB3D_utilities.f90
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
    # X_utilities only needs num_utilities (the index helper c)
    list(APPEND PB3D_CORE_SOURCES Modules/X_utilities.f90)
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
    Modules/flat_vars.f90
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

# include/ and the .mod files of the modules themselves come first; the
# external libraries' include directories come with pb3d_deps
target_include_directories(pb3d_modules PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/mod/full
)
target_link_libraries(pb3d_modules PUBLIC foul pb3d_deps)

#------------------------------------------------------------------------------
# PB3D and POST executables: each is its main program plus the objects of
# pb3d_modules
#------------------------------------------------------------------------------
foreach(prog IN ITEMS PB3D POST)
    add_executable(${prog}
        ${prog}.f90
        $<TARGET_OBJECTS:pb3d_modules>
    )
    target_include_directories(${prog} PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/mod/full
    )
    target_link_libraries(${prog} PRIVATE dfftpack foul pb3d_deps)
endforeach()

#------------------------------------------------------------------------------
# Installation. The installed executables keep the paths of the shared
# libraries they were linked with (PETSc, SLEPc, HDF5, ... from a Spack view
# or a prefix outside the system paths), so they run without LD_LIBRARY_PATH.
#------------------------------------------------------------------------------
set_target_properties(PB3D POST PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
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
    message(STATUS "  PETSc:          ${PETSC_VERSION} (${PETSC_PREFIX})")
    message(STATUS "  distributed LU: ${PB3D_PETSC_DIST_LU}")
    message(STATUS "  SLEPc:          ${SLEPC_VERSION}")
    message(STATUS "  HDF5:           ${HDF5_VERSION} (parallel: ${HDF5_IS_PARALLEL})")
    message(STATUS "  NetCDF-Fortran: ${NETCDF_FORTRAN_VERSION}")
    message(STATUS "  ScaLAPACK:      ${SCALAPACK_FOUND}")
    message(STATUS "  PSPLINE:        ${PSPLINE_LIBRARIES}")
    message(STATUS "  LIBSTELL:       ${LIBSTELL_LIBRARIES}")
    message(STATUS "")
endif()
