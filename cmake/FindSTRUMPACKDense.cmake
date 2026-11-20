# FindSTRUMPACKDense.cmake
# ------------------------
# Find STRUMPACK-Dense 1.1.1, the old dense-only STRUMPACK with a Fortran
# interface. It is not the STRUMPACK 5 to 8 sparse solver, whose interface is
# different. The vacuum module solves its boundary-element system with the
# compressed (HSS) solver of STRUMPACK-Dense when it is found, and with
# ScaLAPACK's LU otherwise.
#
# Looked for only under STRUMPACK_DIR (CMake variable, else environment
# variable, else $HOME/Code/STRUMPACK-Dense-1.1.1) and never in the default
# paths, so that a newer STRUMPACK, from Spack for example, is not picked up.
# The layout of Documentation/spack-setup.md:
#   <STRUMPACK_DIR>/lib or src:      libstrumpack.a
#   <STRUMPACK_DIR>/inc or examples: strumpackdensepackage.mod
#
# Defines:
#   STRUMPACKDense_FOUND    - TRUE if both files were found
#   STRUMPACK_LIBRARY       - the library (cache)
#   STRUMPACK_INCLUDE_DIR   - the directory of the .mod file (cache)
#   STRUMPACK::Dense        - imported target carrying both

if(NOT STRUMPACK_DIR)
    if(DEFINED ENV{STRUMPACK_DIR})
        set(STRUMPACK_DIR "$ENV{STRUMPACK_DIR}")
    else()
        set(STRUMPACK_DIR "$ENV{HOME}/Code/STRUMPACK-Dense-1.1.1")
    endif()
endif()

find_library(STRUMPACK_LIBRARY NAMES strumpack
    PATHS ${STRUMPACK_DIR}/lib ${STRUMPACK_DIR}/src
    NO_DEFAULT_PATH
    DOC "STRUMPACK-Dense 1.1.1 library")
find_path(STRUMPACK_INCLUDE_DIR NAMES strumpackdensepackage.mod
    PATHS ${STRUMPACK_DIR}/inc ${STRUMPACK_DIR}/examples
    NO_DEFAULT_PATH
    DOC "directory of strumpackdensepackage.mod")
mark_as_advanced(STRUMPACK_LIBRARY STRUMPACK_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(STRUMPACKDense
    REQUIRED_VARS STRUMPACK_LIBRARY STRUMPACK_INCLUDE_DIR)

if(STRUMPACKDense_FOUND AND NOT TARGET STRUMPACK::Dense)
    add_library(STRUMPACK::Dense UNKNOWN IMPORTED)
    set_target_properties(STRUMPACK::Dense PROPERTIES
        IMPORTED_LOCATION "${STRUMPACK_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${STRUMPACK_INCLUDE_DIR}")
endif()
