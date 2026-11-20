# PB3DCheckPETSc.cmake
# --------------------
# Checks the PETSc and SLEPc that pkg_check_modules found (the PETSC_* and
# SLEPC_* variables) against what PB3D needs. The configuration of PETSc is
# read from its petscconf.h.
#
#   complex scalars        required: PB3D's matrices are complex
#   double precision       required
#   version 3.19 or newer  required: Modules/SLEPC_ops.f90 handles 3.19 to
#                          the newest tested release with version guards
#   SLEPc minor = PETSc's  required: SLEPc x.y only builds against PETSc x.y
#   newer than tested      warning
#   distributed LU         warning without one: the shift-and-invert
#                          eigensolve factorises a distributed matrix with
#                          LU, and PETSc's own LU works on one process only,
#                          so every PB3D run on two or more MPI processes
#                          stops in EPSSolve (PB3D error 92)
#   64-bit indices         warning: never tested
#
# Sets PB3D_PETSC_DIST_LU to the distributed LU package that PETSc has
# ("MUMPS <version>", or "SuperLU_DIST"), to "none", or to "unknown" when
# petscconf.h was not found.

# the newest PETSc/SLEPc minor release that PB3D was built and tested with
set(PB3D_PETSC_NEWEST_TESTED 3.25)
set(PB3D_PETSC_OLDEST 3.19)

string(REGEX MATCH "^[0-9]+\\.[0-9]+" _pb3d_petsc_mm "${PETSC_VERSION}")
string(REGEX MATCH "^[0-9]+\\.[0-9]+" _pb3d_slepc_mm "${SLEPC_VERSION}")

if(_pb3d_petsc_mm VERSION_LESS PB3D_PETSC_OLDEST)
    message(FATAL_ERROR "PETSc ${PETSC_VERSION} is too old: PB3D needs \
${PB3D_PETSC_OLDEST} or newer.")
endif()
if(NOT _pb3d_slepc_mm VERSION_EQUAL _pb3d_petsc_mm)
    message(FATAL_ERROR "SLEPc ${SLEPC_VERSION} does not belong to PETSc \
${PETSC_VERSION}: their major and minor versions must be equal. Check \
PKG_CONFIG_PATH (it decides which PETSc.pc and SLEPc.pc are found).")
endif()
if(_pb3d_petsc_mm VERSION_GREATER PB3D_PETSC_NEWEST_TESTED)
    message(WARNING "PETSc ${PETSC_VERSION} is newer than the newest release \
PB3D was tested with (${PB3D_PETSC_NEWEST_TESTED}). PETSc changes its Fortran \
interface from time to time: if the build fails in SLEPC_ops.f90, use \
${PB3D_PETSC_NEWEST_TESTED}.")
endif()

# petscconf.h sits in the include directory of the installation, or of
# PETSC_ARCH for a build in place; both are in PETSc.pc
find_file(PB3D_PETSCCONF_H petscconf.h
    PATHS ${PETSC_INCLUDE_DIRS} ${PETSC_INCLUDEDIR}
    NO_DEFAULT_PATH)
mark_as_advanced(PB3D_PETSCCONF_H)

set(PB3D_PETSC_DIST_LU "unknown")
if(NOT PB3D_PETSCCONF_H)
    message(WARNING "petscconf.h not found in the include directories of \
PETSc.pc (${PETSC_INCLUDE_DIRS}): PB3D cannot check that this PETSc has \
complex scalars and a distributed LU.")
else()
    file(STRINGS "${PB3D_PETSCCONF_H}" _pb3d_petscconf
        REGEX "^#define PETSC_(USE_COMPLEX|USE_REAL_DOUBLE|USE_64BIT_INDICES|HAVE_MUMPS|HAVE_SUPERLU_DIST) ")
    macro(_pb3d_petsc_has var macro)
        if(_pb3d_petscconf MATCHES "#define ${macro} 1")
            set(${var} TRUE)
        else()
            set(${var} FALSE)
        endif()
    endmacro()
    _pb3d_petsc_has(_pb3d_complex PETSC_USE_COMPLEX)
    _pb3d_petsc_has(_pb3d_double PETSC_USE_REAL_DOUBLE)
    _pb3d_petsc_has(_pb3d_idx64 PETSC_USE_64BIT_INDICES)
    _pb3d_petsc_has(_pb3d_mumps PETSC_HAVE_MUMPS)
    _pb3d_petsc_has(_pb3d_superlu_dist PETSC_HAVE_SUPERLU_DIST)

    if(NOT _pb3d_complex)
        message(FATAL_ERROR "The PETSc found (${PB3D_PETSCCONF_H}) has real \
scalars; PB3D needs complex ones. Configure PETSc with \
--with-scalar-type=complex (Spack: petsc +complex).")
    endif()
    if(NOT _pb3d_double)
        message(FATAL_ERROR "The PETSc found (${PB3D_PETSCCONF_H}) is not in \
double precision; PB3D needs --with-precision=double.")
    endif()
    if(_pb3d_idx64)
        message(WARNING "This PETSc has 64-bit indices \
(--with-64-bit-indices); PB3D has never been tested with them.")
    endif()
    if(_pb3d_mumps)
        set(PB3D_PETSC_DIST_LU "MUMPS")
        get_filename_component(_pb3d_petsc_inc "${PB3D_PETSCCONF_H}" DIRECTORY)
        set(_pb3d_mumps_ver "")
        if(EXISTS "${_pb3d_petsc_inc}/petscpkg_version.h")
            file(STRINGS "${_pb3d_petsc_inc}/petscpkg_version.h" _pb3d_mumps_ver
                REGEX "^#define PETSC_PKG_MUMPS_VERSION_(MAJOR|MINOR|SUBMINOR) ")
        endif()
        string(REGEX REPLACE ".*MAJOR ([0-9]+).*MINOR ([0-9]+).*SUBMINOR ([0-9]+).*"
            "\\1.\\2.\\3" _pb3d_mumps_ver "${_pb3d_mumps_ver}")
        if(_pb3d_mumps_ver MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
            set(PB3D_PETSC_DIST_LU "MUMPS ${_pb3d_mumps_ver}")
        endif()
    elseif(_pb3d_superlu_dist)
        set(PB3D_PETSC_DIST_LU "SuperLU_DIST")
    else()
        set(PB3D_PETSC_DIST_LU "none")
        message(WARNING "This PETSc has no distributed LU (neither MUMPS nor \
SuperLU_DIST). PB3D builds and runs on one MPI process, but every run on two \
or more processes stops in EPSSolve with error 92. Rebuild PETSc with \
--download-mumps (Spack: petsc +mumps).")
    endif()
endif()
