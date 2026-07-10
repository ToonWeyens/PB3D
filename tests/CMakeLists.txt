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
# Some sandboxed/container environments break OpenMPI's default one-sided
# communication (MPI_Win_create fails in shared-memory/CMA paths, even for
# single-process runs). If this OpenMPI still ships the pt2pt osc component,
# force it for the tests that run the executables; it is slow but always
# works. OpenMPI >= 5 dropped pt2pt (and generally works out of the box), and
# other MPI implementations ignore the variable.
set(PB3D_TEST_MPI_ENV "")
find_program(OMPI_INFO_EXECUTABLE ompi_info)
if(OMPI_INFO_EXECUTABLE)
    execute_process(COMMAND ${OMPI_INFO_EXECUTABLE} --parsable
        OUTPUT_VARIABLE _ompi_info_out ERROR_QUIET)
    if(_ompi_info_out MATCHES "mca:osc:pt2pt")
        set(PB3D_TEST_MPI_ENV "OMPI_MCA_osc=pt2pt")
    endif()
endif()

# Invoking without arguments has to print the usage message. The exit code is
# deliberately not relied upon (PASS_REGULAR_EXPRESSION takes precedence), as
# the error path aborts through MPI.
if(TARGET PB3D)
    add_test(NAME smoke_PB3D_usage COMMAND PB3D)
    set_tests_properties(smoke_PB3D_usage PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage:"
        ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
        LABELS "smoke")
endif()
if(TARGET POST)
    add_test(NAME smoke_POST_usage COMMAND POST)
    set_tests_properties(smoke_POST_usage PROPERTIES
        PASS_REGULAR_EXPRESSION "Usage:"
        ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
        LABELS "smoke")
endif()

#------------------------------------------------------------------------------
# Full-stack tests
#------------------------------------------------------------------------------
# These link against the complete pb3d_modules library and therefore need the
# full dependency stack; they currently cover the vacuum module (Green's
# function kernels against brute-force references, Green's identity on a
# synthetic circular boundary, and the boundary potential solve). Single
# MPI process (the executable runs in MPI singleton mode).
if(PB3D_BUILD_EXECUTABLES)
    add_executable(pb3d_fullstack_tests
        fullstack/main.f90
        fullstack/test_vac_kernels.f90
        fullstack/test_vac_greens.f90
        fullstack/test_vac_3d.f90
        $<TARGET_OBJECTS:pb3d_modules>
    )
    set_target_properties(pb3d_fullstack_tests PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/tests_fullstack)
    target_include_directories(pb3d_fullstack_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/mod/full
    )
    target_link_libraries(pb3d_fullstack_tests PRIVATE
        testdrive
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
            target_link_libraries(pb3d_fullstack_tests PRIVATE PkgConfig::SCALAPACK)
        else()
            target_link_libraries(pb3d_fullstack_tests PRIVATE ${SCALAPACK_LIBRARIES})
        endif()
    endif()
    if(STRUMPACK_FOUND)
        target_link_libraries(pb3d_fullstack_tests PRIVATE ${STRUMPACK_LIBRARIES})
    endif()
    if(LAPACK_FOUND)
        target_link_libraries(pb3d_fullstack_tests PRIVATE ${LAPACK_LIBRARIES})
    endif()
    if(BLAS_FOUND)
        target_link_libraries(pb3d_fullstack_tests PRIVATE ${BLAS_LIBRARIES})
    endif()
    target_link_libraries(pb3d_fullstack_tests PRIVATE m stdc++)

    set(PB3D_FULLSTACK_SUITES
        vac_kernels
        vac_greens
        vac_3d
    )
    foreach(suite IN LISTS PB3D_FULLSTACK_SUITES)
        add_test(NAME fullstack_${suite} COMMAND pb3d_fullstack_tests ${suite})
        set_tests_properties(fullstack_${suite} PROPERTIES
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            LABELS "fullstack")
    endforeach()
endif()

#------------------------------------------------------------------------------
# Physics regression tests (opt-in: need local equilibrium fixtures)
#------------------------------------------------------------------------------
# Pass -DPB3D_FIXTURE_DIR=<dir> at configure time to enable. The fixtures are
# large binary equilibrium files kept outside the repository (see
# Documentation/testing.md for provenance).
if(TARGET PB3D AND PB3D_FIXTURE_DIR)
    if(EXISTS ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12)
        add_test(NAME regression_cbm18a_fixed_bnd
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/cbm18a_fixed_bnd.sh
                $<TARGET_FILE:PB3D>
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_cbm18a_fixed_bnd
                ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12)
        set_tests_properties(regression_cbm18a_fixed_bnd PROPERTIES
            LABELS "regression"
            PASS_REGULAR_EXPRESSION "PASS"
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            TIMEOUT 600)
        add_test(NAME regression_cbm18a_free_bnd
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/cbm18a_free_bnd.sh
                $<TARGET_FILE:PB3D>
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_cbm18a_free_bnd
                ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12)
        set_tests_properties(regression_cbm18a_free_bnd PROPERTIES
            LABELS "regression"
            PASS_REGULAR_EXPRESSION "PASS"
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            TIMEOUT 600)
    else()
        message(STATUS "PB3D_FIXTURE_DIR set but cbm18a_extended_vac.12 not found - regression test skipped")
    endif()
endif()
