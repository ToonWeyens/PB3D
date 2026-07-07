#------------------------------------------------------------------------------
# PB3D test suite
#
# Two layers:
#  - unit tests: build against the dependency-light pb3d_core object library
#    and the vendored test-drive framework; they require no PETSc/SLEPc/HDF5
#    stack and thus run anywhere a Fortran compiler is available.
#  - smoke tests: run the PB3D/POST executables themselves; only registered
#    when those targets are built (PB3D_BUILD_EXECUTABLES=ON).
#
# Run with:  ctest [-L unit|smoke] [-R <test-name>] [--output-on-failure]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Vendored test-drive unit testing framework (see vendor/test-drive/README.md)
#------------------------------------------------------------------------------
add_library(testdrive STATIC vendor/test-drive/testdrive.F90)
set_target_properties(testdrive PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/testdrive)
target_include_directories(testdrive PUBLIC ${CMAKE_BINARY_DIR}/mod/testdrive)

#------------------------------------------------------------------------------
# Unit tests
#------------------------------------------------------------------------------
add_executable(pb3d_unit_tests
    unit/main.f90
    unit/test_str_utilities.f90
    unit/test_files_utilities.f90
    unit/test_dtorh.f90
)
set_target_properties(pb3d_unit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/tests)
target_link_libraries(pb3d_unit_tests PRIVATE pb3d_core testdrive)

# one CTest entry per suite, so failures show up individually and suites can
# be run in parallel
set(PB3D_UNIT_SUITES
    str_utilities
    files_utilities
    dtorh
)
foreach(suite IN LISTS PB3D_UNIT_SUITES)
    add_test(NAME unit_${suite} COMMAND pb3d_unit_tests ${suite})
    set_tests_properties(unit_${suite} PROPERTIES LABELS "unit")
endforeach()

#------------------------------------------------------------------------------
# Smoke tests of the executables
#------------------------------------------------------------------------------
# Invoking without arguments has to print the usage message. The exit code is
# deliberately not relied upon (PASS_REGULAR_EXPRESSION takes precedence), as
# the error path aborts through MPI.
if(TARGET PB3D)
    add_test(NAME smoke_PB3D_usage COMMAND PB3D)
    set_tests_properties(smoke_PB3D_usage PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage:"
        LABELS "smoke")
endif()
if(TARGET POST)
    add_test(NAME smoke_POST_usage COMMAND POST)
    set_tests_properties(smoke_POST_usage PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage:"
        LABELS "smoke")
endif()
