# FindPSPLINE.cmake
# ------------------
# Find the PSPLINE library (Princeton Spline Library)
#
# This module defines:
#   PSPLINE_FOUND        - True if PSPLINE was found
#   PSPLINE_INCLUDE_DIRS - Include directories for PSPLINE
#   PSPLINE_LIBRARIES    - Libraries to link against
#
# The following variables can be set to guide the search:
#   PSPLINE_DIR          - Root directory of PSPLINE installation
#   PSPLINE_ROOT         - Same as PSPLINE_DIR (alternative name)
#   ENV{PSPLINE_DIR}     - Environment variable for PSPLINE root directory

# Use PSPLINE_ROOT or PSPLINE_DIR
set(_PSPLINE_SEARCH_DIRS)
if(PSPLINE_DIR)
    list(APPEND _PSPLINE_SEARCH_DIRS ${PSPLINE_DIR})
endif()
if(PSPLINE_ROOT)
    list(APPEND _PSPLINE_SEARCH_DIRS ${PSPLINE_ROOT})
endif()
if(DEFINED ENV{PSPLINE_DIR})
    list(APPEND _PSPLINE_SEARCH_DIRS $ENV{PSPLINE_DIR})
endif()
if(DEFINED ENV{PSPLINE_ROOT})
    list(APPEND _PSPLINE_SEARCH_DIRS $ENV{PSPLINE_ROOT})
endif()

# Find library
find_library(PSPLINE_LIBRARY
    NAMES pspline
    HINTS ${_PSPLINE_SEARCH_DIRS}
    PATH_SUFFIXES lib lib64 LINUX/lib
    DOC "PSPLINE library"
)

# Find include directory (look for a typical PSPLINE module file)
find_path(PSPLINE_INCLUDE_DIR
    NAMES pspline.mod ezcdf.mod ezspline_obj.mod
    HINTS ${_PSPLINE_SEARCH_DIRS}
    PATH_SUFFIXES include mod LINUX/mod
    DOC "PSPLINE include directory"
)

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PSPLINE
    REQUIRED_VARS PSPLINE_LIBRARY PSPLINE_INCLUDE_DIR
    FAIL_MESSAGE "Could not find PSPLINE. Set PSPLINE_DIR to the installation root."
)

# Set output variables
if(PSPLINE_FOUND)
    set(PSPLINE_LIBRARIES ${PSPLINE_LIBRARY})
    set(PSPLINE_INCLUDE_DIRS ${PSPLINE_INCLUDE_DIR})

    # Create imported target
    if(NOT TARGET PSPLINE::PSPLINE)
        add_library(PSPLINE::PSPLINE UNKNOWN IMPORTED)
        set_target_properties(PSPLINE::PSPLINE PROPERTIES
            IMPORTED_LOCATION "${PSPLINE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PSPLINE_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(PSPLINE_LIBRARY PSPLINE_INCLUDE_DIR)
