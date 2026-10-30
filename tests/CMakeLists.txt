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
set(PB3D_UNIT_SOURCES
    unit/main.f90
    unit/test_str_utilities.f90
    unit/test_files_utilities.f90
    unit/test_dtorh.f90
    unit/test_PB3D_utilities.f90
)
set(PB3D_UNIT_SUITES
    str_utilities
    files_utilities
    dtorh
    PB3D_utilities
)
# the numerics suites need num_utilities/num_ops in pb3d_core (see the LAPACK
# and ldebug conditions in the top-level CMakeLists.txt); X_utilities uses
# num_utilities's index helper c
if(PB3D_CORE_HAS_NUM)
    list(APPEND PB3D_UNIT_SOURCES unit/test_num_utilities.f90
        unit/test_X_utilities.f90)
    list(APPEND PB3D_UNIT_SUITES num_utilities X_utilities)
endif()
if(PB3D_CORE_HAS_NUM_OPS)
    list(APPEND PB3D_UNIT_SOURCES unit/test_num_ops.f90)
    list(APPEND PB3D_UNIT_SUITES num_ops)
endif()

add_executable(pb3d_unit_tests ${PB3D_UNIT_SOURCES})
set_target_properties(pb3d_unit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/tests)
target_link_libraries(pb3d_unit_tests PRIVATE pb3d_core testdrive)
if(PB3D_CORE_HAS_NUM)
    target_compile_definitions(pb3d_unit_tests PRIVATE PB3D_TEST_NUM)
endif()
if(PB3D_CORE_HAS_NUM_OPS)
    target_compile_definitions(pb3d_unit_tests PRIVATE PB3D_TEST_NUM_OPS)
endif()

# one CTest entry per suite, so failures show up individually and suites can
# be run in parallel
foreach(suite IN LISTS PB3D_UNIT_SUITES)
    add_test(NAME unit_${suite} COMMAND pb3d_unit_tests ${suite})
    set_tests_properties(unit_${suite} PROPERTIES LABELS "unit")
endforeach()

# index logic of setup_modes (X_ops) replayed in python on the mode-window
# sequences that broke the slot bookkeeping (issue #4, step 5)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME unit_setup_modes_sim
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reference/setup_modes_sim.py)
    set_tests_properties(unit_setup_modes_sim PROPERTIES
        PASS_REGULAR_EXPRESSION "PASS"
        LABELS "unit")
endif()

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
# synthetic circular boundary, and the boundary potential solve). The suites
# are rank-agnostic (they assert on globally gathered quantities, see
# fullstack/fullstack_utils.f90): they run in MPI singleton mode by default
# and are additionally registered under mpirun -n 2 and 4 below, which
# exercises the 2-D block-cyclic ScaLAPACK distribution (blocksize 16)
# instead of the single whole-matrix block.
if(PB3D_BUILD_EXECUTABLES)
    add_executable(pb3d_fullstack_tests
        fullstack/main.f90
        fullstack/fullstack_utils.f90
        fullstack/test_vac_kernels.f90
        fullstack/test_vac_greens.f90
        fullstack/test_vac_3d.f90
        fullstack/test_splines.f90
        fullstack/test_calc_int_vol.f90
        fullstack/test_read_HEL.f90
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
        splines
        calc_int_vol
    )
    foreach(suite IN LISTS PB3D_FULLSTACK_SUITES)
        add_test(NAME fullstack_${suite} COMMAND pb3d_fullstack_tests ${suite})
        set_tests_properties(fullstack_${suite} PROPERTIES
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            LABELS "fullstack")
    endforeach()

    # Multi-rank rows for the suites with distributed linear algebra
    # (vac_kernels is purely local kernel math, so multi-rank runs of it
    # would only duplicate work). OMPI_ALLOW_RUN_AS_ROOT covers container
    # use and is ignored by non-root/non-OpenMPI setups; --oversubscribe
    # keeps -n 4 working on smaller CI runners.
    set(PB3D_FULLSTACK_MPI_SUITES
        vac_greens
        vac_3d
    )
    if(MPIEXEC_EXECUTABLE)
        set(_mpi_test_env "${PB3D_TEST_MPI_ENV}")
        list(APPEND _mpi_test_env
            "OMPI_ALLOW_RUN_AS_ROOT=1" "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
        set(_mpiexec_flags "")
        if(OMPI_INFO_EXECUTABLE)
            set(_mpiexec_flags "--oversubscribe")
        endif()
        foreach(np IN ITEMS 2 4)
            foreach(suite IN LISTS PB3D_FULLSTACK_MPI_SUITES)
                add_test(NAME fullstack_${suite}_np${np}
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                        ${_mpiexec_flags}
                        $<TARGET_FILE:pb3d_fullstack_tests> ${suite})
                set_tests_properties(fullstack_${suite}_np${np} PROPERTIES
                    ENVIRONMENT "${_mpi_test_env}"
                    PROCESSORS ${np}
                    LABELS "fullstack;mpi")
            endforeach()
        endforeach()
    endif()
endif()

option(PB3D_TEST_MPI_REGRESSION
    "Also run the cbm18a fixed-boundary regression deck under mpirun -n 2 and 4 and require the eigenvalue to agree with the one-process run to 1e-6 (needs a distributed LU in PETSc, e.g. MUMPS)"
    OFF)

#------------------------------------------------------------------------------
# Physics regression tests
#------------------------------------------------------------------------------
# The equilibrium fixtures are committed xz-compressed under tests/fixtures
# (cbm18a_extended_vac.12: 26 MB -> 5 MB) and decompressed into the build
# tree at configure time, so the regression layer runs anywhere - including
# CI - without external files. Passing -DPB3D_FIXTURE_DIR=<dir> overrides
# this with a local fixture directory (see Documentation/testing.md).
if(NOT PB3D_FIXTURE_DIR AND
   EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cbm18a_extended_vac.12.xz)
    find_program(XZ_EXECUTABLE xz)
    if(XZ_EXECUTABLE)
        set(PB3D_FIXTURE_DIR ${CMAKE_BINARY_DIR}/fixtures)
        file(MAKE_DIRECTORY ${PB3D_FIXTURE_DIR})
        # every committed *.12.xz (cbm18a and, since issue #4 step 0, the two
        # circular external-kink cases) and *.nc.xz (the VMEC twin of the
        # eps = 0.33 circular case, issue #4 step 6e)
        file(GLOB _fixtures_xz ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/*.12.xz
                               ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/*.nc.xz)
        foreach(_fx ${_fixtures_xz})
            get_filename_component(_fx_name ${_fx} NAME_WLE)
            if(NOT EXISTS ${PB3D_FIXTURE_DIR}/${_fx_name})
                message(STATUS "Decompressing regression fixture ${_fx_name} to ${PB3D_FIXTURE_DIR}")
                execute_process(
                    COMMAND ${XZ_EXECUTABLE} -dkc ${_fx}
                    OUTPUT_FILE ${PB3D_FIXTURE_DIR}/${_fx_name}
                    RESULT_VARIABLE _xz_res)
                if(NOT _xz_res EQUAL 0)
                    message(WARNING "Fixture decompression failed - regression tests disabled")
                    unset(PB3D_FIXTURE_DIR)
                    break()
                endif()
            endif()
        endforeach()
    else()
        message(STATUS "xz not found - regression tests disabled (or pass PB3D_FIXTURE_DIR)")
    endif()
endif()

# the read_HEL golden-file suite needs the fixture, so it is registered
# here rather than with the other fullstack suites above
if(TARGET pb3d_fullstack_tests AND PB3D_FIXTURE_DIR)
    if(EXISTS ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12)
        add_test(NAME fullstack_read_HEL
            COMMAND pb3d_fullstack_tests read_HEL)
        set_tests_properties(fullstack_read_HEL PROPERTIES
            ENVIRONMENT "${PB3D_TEST_MPI_ENV};PB3D_FIXTURE_DIR=${PB3D_FIXTURE_DIR}"
            LABELS "fullstack")
    endif()
endif()

if(TARGET PB3D AND PB3D_FIXTURE_DIR)
    # the HDF5 checks (shear, vacuum response matrix) need a python3 with h5py;
    # pass -DPython3_EXECUTABLE=... to pick one. Without it the eigenvalue
    # assertions still run and the HDF5 checks are skipped.
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import h5py"
            RESULT_VARIABLE _h5py_res OUTPUT_QUIET ERROR_QUIET)
    else()
        set(_h5py_res 1)
    endif()
    if(_h5py_res EQUAL 0)
        set(PB3D_TEST_PYTHON ${Python3_EXECUTABLE})
    else()
        set(PB3D_TEST_PYTHON "")
        message(STATUS "python3 with h5py not found - the HDF5 regression checks are skipped")
    endif()
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
                ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12
                ${PB3D_TEST_PYTHON})
        set_tests_properties(regression_cbm18a_free_bnd PROPERTIES
            LABELS "regression"
            PASS_REGULAR_EXPRESSION "PASS"
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            TIMEOUT 600)
        # the shear check reads PB3D_out.h5 with h5py
        if(PB3D_TEST_PYTHON)
            add_test(NAME regression_cbm18a_shear_HEL
                COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/cbm18a_shear_HEL.sh
                    $<TARGET_FILE:PB3D>
                    ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_cbm18a_fixed_bnd
                    ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12
                    ${Python3_EXECUTABLE})
            set_tests_properties(regression_cbm18a_shear_HEL PROPERTIES
                LABELS "regression"
                PASS_REGULAR_EXPRESSION "PASS"
                ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
                TIMEOUT 600)
        endif()
        # multi-process consistency of the eigenvalue (issue #4: the 2e-5 shift of
        # redistribute_var sat inside the 1e-3 anchor tolerance); opt-in because
        # more than one process needs a distributed LU inside PETSc
        if(PB3D_TEST_MPI_REGRESSION AND MPIEXEC_EXECUTABLE)
            foreach(np IN ITEMS 2 4)
                add_test(NAME regression_cbm18a_fixed_bnd_np${np}
                    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/cbm18a_fixed_bnd_np.sh
                        $<TARGET_FILE:PB3D>
                        ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_cbm18a_fixed_bnd
                        ${PB3D_FIXTURE_DIR}/cbm18a_extended_vac.12
                        "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${_mpiexec_flags}")
                set_tests_properties(regression_cbm18a_fixed_bnd_np${np} PROPERTIES
                    LABELS "regression;mpi"
                    PASS_REGULAR_EXPRESSION "PASS"
                    ENVIRONMENT "${_mpi_test_env}"
                    PROCESSORS ${np}
                    TIMEOUT 900)
            endforeach()
        endif()
    else()
        message(STATUS "PB3D_FIXTURE_DIR set but cbm18a_extended_vac.12 not found - regression test skipped")
    endif()
    # the circular external kink against MISHKA (issue #4, steps 0 and 3)
    if(EXISTS ${PB3D_FIXTURE_DIR}/circ_eps005_q104.12)
        add_test(NAME regression_circ_eps005_free_bnd
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/circ_eps005_free_bnd.sh
                $<TARGET_FILE:PB3D>
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_circ_eps005
                ${PB3D_FIXTURE_DIR}/circ_eps005_q104.12
                ${PB3D_TEST_PYTHON})
        set_tests_properties(regression_circ_eps005_free_bnd PROPERTIES
            LABELS "regression"
            PASS_REGULAR_EXPRESSION "PASS"
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            TIMEOUT 600)
    endif()
    # the HELENA fixture against its VMEC twin, fixed boundary (issue #4, step 6e)
    if(EXISTS ${PB3D_FIXTURE_DIR}/circ_eps033_q104.12 AND
       EXISTS ${PB3D_FIXTURE_DIR}/wout_circ_eps033_q104.nc)
        add_test(NAME regression_circ_eps033_vmec_twin
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/regression/circ_eps033_vmec_twin.sh
                $<TARGET_FILE:PB3D>
                ${CMAKE_CURRENT_SOURCE_DIR}/regression/input_circ_eps033
                ${PB3D_FIXTURE_DIR}/circ_eps033_q104.12
                ${PB3D_FIXTURE_DIR}/wout_circ_eps033_q104.nc)
        set_tests_properties(regression_circ_eps033_vmec_twin PROPERTIES
            LABELS "regression"
            PASS_REGULAR_EXPRESSION "PASS"
            ENVIRONMENT "${PB3D_TEST_MPI_ENV}"
            TIMEOUT 600)
    endif()
endif()
