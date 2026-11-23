# FindSCALAPACK.cmake
# -------------------
# Find ScaLAPACK together with its BLACS. PB3D calls both directly: the BLACS
# process grid of the vacuum module (vac_vars.f90, MPI_ops.f90) and the LU
# solve of the vacuum boundary-element system. ScaLAPACK is therefore
# required.
#
# Search order:
#   1. SCALAPACK_LIBRARIES given on the command line: used as it is (for
#      example the MKL pair mkl_scalapack_lp64;mkl_blacs_openmpi_lp64).
#   2. pkg-config, under the names scalapack (Spack, the Snellius module),
#      scalapack-openmpi and scalapack-mpich (Debian, Ubuntu).
#   3. find_library under SCALAPACK_DIR or SCALAPACK_ROOT (CMake or
#      environment variable), the library directories of PETSc (a PETSc
#      configured with --download-scalapack installs it there), and the
#      default paths.
#
# Defines:
#   SCALAPACK_FOUND        - TRUE if ScaLAPACK was found
#   SCALAPACK_LIBRARIES    - what to link
#   SCALAPACK::SCALAPACK   - imported target carrying SCALAPACK_LIBRARIES

include(FindPackageHandleStandardArgs)

if(SCALAPACK_LIBRARIES)
    set(_SCALAPACK_HOW "SCALAPACK_LIBRARIES given")
else()
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_search_module(_SCALAPACK_PC QUIET IMPORTED_TARGET
            scalapack scalapack-openmpi scalapack-mpich)
    endif()
    if(_SCALAPACK_PC_FOUND)
        set(SCALAPACK_LIBRARIES PkgConfig::_SCALAPACK_PC)
        set(_SCALAPACK_HOW "pkg-config ${_SCALAPACK_PC_MODULE_NAME}")
    else()
        find_library(SCALAPACK_LIBRARY
            NAMES scalapack scalapack-openmpi scalapack-mpich
            HINTS ${SCALAPACK_DIR} ${SCALAPACK_ROOT}
                  $ENV{SCALAPACK_DIR} $ENV{SCALAPACK_ROOT}
                  ${PETSC_LIBRARY_DIRS}
            PATH_SUFFIXES lib lib64
            DOC "ScaLAPACK library")
        mark_as_advanced(SCALAPACK_LIBRARY)
        if(SCALAPACK_LIBRARY)
            set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
            set(_SCALAPACK_HOW "${SCALAPACK_LIBRARY}")
        endif()
    endif()
endif()

find_package_handle_standard_args(SCALAPACK
    REQUIRED_VARS SCALAPACK_LIBRARIES
    REASON_FAILURE_MESSAGE "PB3D needs ScaLAPACK (with BLACS). Install it \
(Spack: netlib-scalapack, Ubuntu: libscalapack-openmpi-dev), or point \
SCALAPACK_DIR at its prefix, or give the libraries in SCALAPACK_LIBRARIES.")

if(SCALAPACK_FOUND AND NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED)
    set_target_properties(SCALAPACK::SCALAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES "${SCALAPACK_LIBRARIES}")
    message(STATUS "ScaLAPACK: ${_SCALAPACK_HOW}")
endif()
