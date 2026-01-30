# FindLIBSTELL.cmake
# ------------------
# Find the LIBSTELL library from STELLOPT
#
# This module defines:
#   LIBSTELL_FOUND        - True if LIBSTELL was found
#   LIBSTELL_INCLUDE_DIRS - Include directories for LIBSTELL
#   LIBSTELL_LIBRARIES    - Libraries to link against
#
# The following variables can be set to guide the search:
#   LIBSTELL_DIR          - Root directory of LIBSTELL installation
#   LIBSTELL_ROOT         - Same as LIBSTELL_DIR (alternative name)
#   ENV{LIBSTELL_DIR}     - Environment variable for LIBSTELL root directory
#   STELLOPT_DIR          - Root directory of STELLOPT (LIBSTELL is usually in bin/)

# Use LIBSTELL_ROOT, LIBSTELL_DIR, or STELLOPT_DIR
set(_LIBSTELL_SEARCH_DIRS)
if(LIBSTELL_DIR)
    list(APPEND _LIBSTELL_SEARCH_DIRS ${LIBSTELL_DIR})
endif()
if(LIBSTELL_ROOT)
    list(APPEND _LIBSTELL_SEARCH_DIRS ${LIBSTELL_ROOT})
endif()
if(STELLOPT_DIR)
    list(APPEND _LIBSTELL_SEARCH_DIRS ${STELLOPT_DIR})
    list(APPEND _LIBSTELL_SEARCH_DIRS ${STELLOPT_DIR}/bin)
endif()
if(DEFINED ENV{LIBSTELL_DIR})
    list(APPEND _LIBSTELL_SEARCH_DIRS $ENV{LIBSTELL_DIR})
endif()
if(DEFINED ENV{LIBSTELL_ROOT})
    list(APPEND _LIBSTELL_SEARCH_DIRS $ENV{LIBSTELL_ROOT})
endif()
if(DEFINED ENV{STELLOPT_DIR})
    list(APPEND _LIBSTELL_SEARCH_DIRS $ENV{STELLOPT_DIR})
    list(APPEND _LIBSTELL_SEARCH_DIRS $ENV{STELLOPT_DIR}/bin)
endif()

# Find library
find_library(LIBSTELL_LIBRARY
    NAMES stell libstell
    HINTS ${_LIBSTELL_SEARCH_DIRS}
    PATH_SUFFIXES lib lib64 bin
    DOC "LIBSTELL library"
)

# Find include directory (look for libstell module directory)
# LIBSTELL puts its module files in a libstell_dir subdirectory
find_path(LIBSTELL_INCLUDE_DIR
    NAMES vsvd0.mod vmec_input.mod
    HINTS ${_LIBSTELL_SEARCH_DIRS}
    PATH_SUFFIXES include mod libstell_dir
    DOC "LIBSTELL include directory"
)

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBSTELL
    REQUIRED_VARS LIBSTELL_LIBRARY LIBSTELL_INCLUDE_DIR
    FAIL_MESSAGE "Could not find LIBSTELL. Set LIBSTELL_DIR to the installation root (or STELLOPT_DIR for STELLOPT)."
)

# Set output variables
if(LIBSTELL_FOUND)
    set(LIBSTELL_LIBRARIES ${LIBSTELL_LIBRARY})
    set(LIBSTELL_INCLUDE_DIRS ${LIBSTELL_INCLUDE_DIR})

    # Create imported target
    if(NOT TARGET LIBSTELL::LIBSTELL)
        add_library(LIBSTELL::LIBSTELL UNKNOWN IMPORTED)
        set_target_properties(LIBSTELL::LIBSTELL PROPERTIES
            IMPORTED_LOCATION "${LIBSTELL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${LIBSTELL_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(LIBSTELL_LIBRARY LIBSTELL_INCLUDE_DIR)
