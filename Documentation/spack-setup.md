# Installing PB3D

This guide installs PB3D on macOS or Linux: the libraries from a Spack
environment (`spack.yaml` in the project root), PSPLINE and LIBSTELL with
the scripts in `Libraries/`, then PB3D itself with CMake. For the Snellius
supercomputer, which uses its own modules and a PETSc built from source, see
`snellius-setup.md`; for the tests, `testing.md`.

What PB3D needs, and what CMake checks at configure time:

| Library | Needed for | Checked |
|---|---|---|
| MPI (Fortran and C) | everything; PB3D is an MPI program | found by `find_package(MPI)` |
| PETSc, complex scalars, double precision, 3.19 or newer | the eigenvalue problem | stops if real, single precision or older than 3.19; warns if newer than 3.25, the newest tested |
| a distributed LU in PETSc (MUMPS or SuperLU_DIST) | runs on two or more processes | warns if absent: PB3D then runs on one process only |
| SLEPc of the same minor version as PETSc | the eigensolver | stops if the versions differ |
| HDF5 with Fortran, built with MPI | the output | stops if serial |
| NetCDF-Fortran | VMEC equilibria, through LIBSTELL | found by pkg-config |
| ScaLAPACK with BLACS | the vacuum | stops if absent |
| PSPLINE | splines | `cmake/FindPSPLINE.cmake` |
| LIBSTELL (`read_wout_mod`) | VMEC equilibria | `cmake/FindLIBSTELL.cmake` |

## 1. Install Spack

```bash
git clone https://github.com/spack/spack.git $HOME/spack
# add to your shell (zsh shown)
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.zshrc
source ~/.zshrc
```

Spack ships its own Python, so no extra dependencies are needed. Verify the CLI:

```bash
spack --version
```

## 2. Bootstrap compilers (macOS)

Apple Clang does not provide `gfortran`, so you need GCC with Fortran support. On macOS, you can install it via Homebrew:

```bash
brew install gcc
```

Then register the compiler with Spack:

```bash
spack compiler find
spack compilers  # Should show gcc with Fortran support
```

## 3. Create a PB3D environment

```bash
# reuse an existing MPI (e.g. Homebrew open-mpi) instead of building one:
spack external find cmake openmpi gmake perl python

cd /path/to/PB3D
spack env create pb3d spack.yaml
spack env activate pb3d
```

The `spack.yaml` in the project root lists the libraries: PETSc 3.19 to
3.25 with complex scalars and MUMPS, the SLEPc that belongs to it, parallel
HDF5, NetCDF-Fortran and ScaLAPACK. Spack takes the newest PETSc in that
range that it knows. MUMPS gives PETSc the distributed LU that PB3D needs
on two or more processes. MUMPS is compiled with `-DAVOID_MPI_IN_PLACE`.
Without it, PB3D on two or more processes crashed on macOS with a bus error
in `zmumps_scatter_dist_rhs`: each shared MUMPS library keeps its own copy
of the COMMON block behind `MPI_IN_PLACE` in `mpif.h`, Open MPI does not
recognise the `MPI_IN_PLACE` that MUMPS passes, and a reduction in MUMPS
adds up garbage. For a PETSc whose MUMPS was built without the flag, the
PETSc option `-st_mat_mumps_icntl_20 0` avoids that reduction: MUMPS then
gathers the right-hand side on one process. OpenMPI is preferred, not
required: on a cluster, `spack external find` the site MPI first and Spack
uses it.

Now concretize and install:

```bash
spack concretize -f
spack install
```

The first build can take some time because PETSc/SLEPc compile a large stack. Spack caches results, so later rebuilds are much faster.
An environment made from an older `spack.yaml` (PETSc without MUMPS) keeps
working on one process; `spack concretize -f && spack install` in it after
updating its `spack.yaml` adds MUMPS and rebuilds PETSc and SLEPc only.

### macOS: `mpirun` of Homebrew's OpenMPI 5

On Apple silicon, `mpirun` of OpenMPI 5 can segfault before any program
starts, in `prte_hwloc_base_get_topo_signature` (seen with OpenMPI 5.0.9 and
hwloc 2.13 on an M3 Pro under macOS 27, whose six efficiency and five
performance cores sit under two L2 caches of different sizes, with no L3).
A synthetic hwloc topology with the same number of cores avoids it. macOS
cannot bind processes to cores anyway, so nothing is lost. Put this in
`~/.zshrc`:

```bash
export HWLOC_SYNTHETIC="package:1 core:$(sysctl -n hw.physicalcpu) pu:1"
```

The CMake configuration of the tests checks whether `mpiexec` needs the
variable and sets it for the tests itself. A run without `mpirun`
(`./PB3D input eq`, MPI singleton mode) never needed it.

## 4. Build PSPLINE (not in Spack)

PSPLINE is a Princeton spline library required by PB3D. Its legacy NTCC
make system does not work with modern toolchains; use the provided script
(needs only gfortran, no netcdf):

```bash
git clone https://github.com/ilonster/pspline.git ~/Code/pspline
cp /path/to/PB3D/Libraries/build_pspline.sh ~/Code/pspline/
~/Code/pspline/build_pspline.sh
export PSPLINE_DIR=$HOME/Code/pspline/install
```

## 5. Build LIBSTELL (not in Spack)

PB3D only uses `read_wout_mod` from LIBSTELL. The recommended route is the
minimal build script, which compiles just that module and its transitive
dependencies from a STELLOPT checkout (needs mpif90 + netcdf-fortran from
the Spack env):

```bash
git clone --depth 1 https://github.com/PrincetonUniversity/STELLOPT.git ~/Code/stellopt
/path/to/PB3D/Libraries/build_libstell_min.sh
export LIBSTELL_DIR=/path/to/PB3D/Libraries/libstell-min
```

Alternatively, build the full LIBSTELL via STELLOPT's own build system:

```bash
git clone https://github.com/PrincetonUniversity/STELLOPT.git
cd STELLOPT

# Create spack-compatible make.inc (macOS example)
cat > SHARE/make_spack.inc << 'EOF'
STELLOPT_HOME ?= $(HOME)/Code/stellopt/bin
MYHOME = $(STELLOPT_HOME)
SHELL = /bin/sh
PWD1 = `pwd`
PRECOMP:= cpp -traditional-cpp -E -P -C -DMACOSX
COMPILE = gfortran
COMPILE_FREE = gfortran -ffree-form -ffree-line-length-none -ffixed-line-length-none
LINK    = ld $(FLAGS) -o
LINK_AR = ar -ruvs
LINK_C  = gcc -shared -Wl,-export-dynamic

# Point to Spack view
SPACK_VIEW = $(shell spack env location --view)
SCALAPACK_HOME = $(SPACK_VIEW)
HDF5_HOME = $(SPACK_VIEW)
NETCDF_HOME = $(SPACK_VIEW)

FLAGS_R = -O2 -g -fexternal-blas -fallow-argument-mismatch
FLAGS_D = -O0 -g -fexternal-blas -fbacktrace -fcheck=all,no-array-temps -fbounds-check -fallow-argument-mismatch
LIBS    = -L$(SPACK_VIEW)/lib -lscalapack -lvecLibFort -framework Accelerate

LMPI    = T
MPI_COMPILE = mpif90
MPI_COMPILE_FREE = mpif90 -ffree-form -ffree-line-length-none -ffixed-line-length-none
MPI_COMPILE_C = mpicc
MPI_LINK = mpif90 -Wl,-no_compact_unwind
MPI_RUN = mpiexec
MPI_RUN_OPTS = --use-hwthread-cpus

LNAG = F
NAG_LIB =

LNETCDF = T
NETCDF_INC = -I$(NETCDF_HOME)/include
NETCDF_LIB = -L$(NETCDF_HOME)/lib -lnetcdff -lnetcdf

LFFTW3 = F
FFTW3_INC =
FFTW3_LIB =

LHDF5 = T
HDF5_INC = -I$(HDF5_HOME)/include
HDF5_LIB = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm

LPGPLOT = F
LSILO = F
LDKES = T
LNEO  = T
LGENE = F
LCOILOPT = F
LTERPSICHORE = F
LTRAVIS = F
LREGCOIL = F
LSFINCS = F
LAEOPT = F
LMANGO = F

LIB_SHARE = -lc -lgfortran -lstdc++ -lmpi -lmpi_mpifh -lz -lc -lm -lpthread $(LIBS) -lc
EOF

# Set MACHINE to use the spack config
export MACHINE=spack
export STELLOPT_PATH=$PWD

# Build LIBSTELL
cd LIBSTELL && make release && cd ..

# The library ends up at $STELLOPT_PATH/install/lib/libstell.a
# Module files are at $STELLOPT_PATH/install/include/*.mod
```

## 6. Building PB3D with CMake

With the Spack environment active and PSPLINE/LIBSTELL built:

```bash
spack env activate pb3d

cmake -S . -B build \
    -DPSPLINE_DIR=$HOME/Code/pspline/install \
    -DLIBSTELL_DIR=$PWD/Libraries/libstell-min \
    -DPB3D_FIXTURE_DIR=$HOME/Code/PB3D-fixtures    # optional, see testing.md
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The executables `PB3D` and `POST` are created in the build directory. The
summary at the end of the configure step names the PETSc it found, its
distributed LU and whether HDF5 is parallel; read it once. `cmake --install build --prefix <dir>` copies the two
executables to `<dir>/bin`; they keep the paths of the Spack libraries, so
they run without `LD_LIBRARY_PATH` or `DYLD_LIBRARY_PATH`.

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `PB3D_ENABLE_DEBUG` | OFF | Enable debug mode (ldebug preprocessor flag) |
| `PB3D_ENABLE_INFINIBAND` | OFF | Enable InfiniBand support (lIB flag) |
| `PB3D_BUILD_EXECUTABLES` | ON | Build PB3D/POST (needs the full stack); OFF builds only the core library + unit tests |
| `BUILD_TESTING` | ON | Build the CTest test suite (see `testing.md`) |
| `CMAKE_BUILD_TYPE` | Release | Build type (Debug, Release, RelWithDebInfo) |
| `PSPLINE_DIR` | - | Path to PSPLINE installation |
| `LIBSTELL_DIR` | - | Path to LIBSTELL installation |
| `SCALAPACK_DIR`, `SCALAPACK_LIBRARIES` | - | Where ScaLAPACK is, when pkg-config does not know it; or the libraries themselves (MKL) |
| `PB3D_FIXTURE_DIR` | build/fixtures | Equilibrium files of the regression tests (`testing.md`) |
| `PB3D_TEST_MPI_REGRESSION` | ON when PETSc has a distributed LU | The 2- and 4-process regression rows |
| `Python3_EXECUTABLE` | `.venv/bin/python` if it exists | Python with h5py, numpy and scipy for the HDF5 checks (`uv sync` makes `.venv`) |

### CMake versions

CMake 3.18 or newer is enough. The Makefile generator of CMake 3.31 used to
hang in its dependency scan of PETSc 3.25's Fortran headers, which include
each other; `CMakeLists.txt` now limits that scan to PB3D's own files
(`CMAKE_DEPENDS_IN_PROJECT_ONLY`), and CMake 3.31 builds normally. After
reinstalling PETSc in place, clean the build directory, because a change in
an external header no longer triggers a rebuild by itself.

### Testing without the full stack

The unit tests need no external libraries at all (see `testing.md`):

```bash
cmake -S . -B build-tests -DPB3D_BUILD_EXECUTABLES=OFF
cmake --build build-tests -j
ctest --test-dir build-tests --output-on-failure
```

### Debug Build

```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPB3D_ENABLE_DEBUG=ON \
    -DPSPLINE_DIR=$PSPLINE_DIR \
    -DLIBSTELL_DIR=$LIBSTELL_DIR
make -j4
```

## 7. HPC Cluster with Module System

On a cluster, the site modules usually give MPI, HDF5, NetCDF-Fortran and
ScaLAPACK. The site PETSc is often built with real scalars, which PB3D
cannot use (CMake stops with a message saying so); then build PETSc with
complex scalars and MUMPS yourself, as `snellius-setup.md` and
`scripts/snellius/build_deps.sh` do on Snellius, or use Spack with the site
MPI registered as external.

```bash
# Load the modules (names differ per site)
module load gcc openmpi hdf5 netcdf-fortran scalapack pkgconf
export PKG_CONFIG_PATH=<petsc-prefix>/lib/pkgconfig:<slepc-prefix>/lib/pkgconfig:$PKG_CONFIG_PATH

cmake -S . -B build \
    -DPSPLINE_DIR=/path/to/pspline \
    -DLIBSTELL_DIR=/path/to/libstell
cmake --build build -j
```

## 8. Legacy Makefile

The original Makefile is preserved as `Makefile.legacy` for reference. To use it:

1. Copy it back: `cp Makefile.legacy Makefile`
2. Edit the paths at the top for your system
3. Run `make all`

The CMake-based build is recommended for new installations.

## 9. Sharing the setup

`spack concretize` writes `spack.lock` next to the environment's
`spack.yaml` (`spack env location`). It pins every version, variant and
hash; someone else reproduces the same stack with

```bash
spack env create pb3d /path/to/spack.lock
```

## 10. Useful Spack commands

- `spack find` – list installed packages.
- `spack find --deps <spec>` – show dependency tree.
- `spack uninstall <spec>` – remove a package.
- `spack clean -a` – clear build caches if disk space is tight.
- `spack env location --view` – get path to environment view

For more background, see the official docs: <https://spack.readthedocs.io/>.
