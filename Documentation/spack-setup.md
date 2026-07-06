# PB3D Spack Environment

This guide walks you through creating a reproducible Spack environment for PB3D on macOS or Linux. It supplies the PETSc/SLEPc/HDF5 toolchain PB3D expects, while keeping all binaries isolated from your host system.

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
cd /path/to/PB3D
spack env create pb3d-env spack.yaml
spack env activate pb3d-env
```

The `spack.yaml` in the project root contains the dependency list. You can edit it before `concretize` if your cluster needs different MPI or math libs.

Now concretize and install:

```bash
spack concretize -f
spack install
```

The first build can take some time because PETSc/SLEPc compile a large stack. Spack caches results, so later rebuilds are much faster.

## 4. Build PSPLINE (not in Spack)

PSPLINE is a Princeton spline library required by PB3D. Download and build it:

```bash
git clone https://github.com/ilonster/pspline.git
cd pspline
# Build with the Spack environment active
export PSPLINE_DIR=$HOME/opt/pspline
mkdir -p $PSPLINE_DIR
# Follow PSPLINE's build instructions, typically:
# make FORTRAN=gfortran
# cp -r lib mod $PSPLINE_DIR/
```

## 5. Build LIBSTELL (not in Spack)

LIBSTELL is part of the STELLOPT suite. Download and build it:

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

## 6. Build STRUMPACK-Dense 1.1.1 (required for vacuum module)

PB3D's vacuum module requires the old STRUMPACK-Dense 1.1.1 library (not the newer STRUMPACK 7.x):

```bash
cd ~/Code
curl -L -O http://portal.nersc.gov/project/sparse/strumpack/STRUMPACK-Dense-1.1.1.tar.gz
tar -xzf STRUMPACK-Dense-1.1.1.tar.gz
cd STRUMPACK-Dense-1.1.1/examples

# Create Makefile.inc for macOS with Spack
cat > Makefile.inc << 'EOF'
SPACK_VIEW = $(shell spack env location --view)
CXXFLAGS  = -O3
CFLAGS    = -O3
FFLAGS    = -O3 -fallow-argument-mismatch
CXX       = mpic++
CC        = mpicc
FC        = mpif90
LIB       = -L$(SPACK_VIEW)/lib -lscalapack -lvecLibFort -framework Accelerate -lm
LIBCXX    = -lstdc++ -lc++
EOF

# Build Fortran interface
make f90_example

# Create library structure
cd ..
mkdir -p lib inc
cd src && ar -rcs libstrumpack.a StrumpackDensePackage_C.o StrumpackDensePackage.o && cd ..
ln -sf ../src/libstrumpack.a lib/libstrumpack.a
ln -sf ../examples/strumpackdensepackage.mod inc/strumpackdensepackage.mod

export STRUMPACK_DIR=$PWD
```

## 7. Building PB3D with CMake

With the Spack environment active and PSPLINE/LIBSTELL/STRUMPACK-Dense built:

```bash
spack env activate pb3d-env

# Set paths to external dependencies
export PSPLINE_DIR=$HOME/Code/pspline/install
export LIBSTELL_DIR=$HOME/Code/stellopt/install
export STRUMPACK_DIR=$HOME/Code/STRUMPACK-Dense-1.1.1

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. \
    -DPSPLINE_DIR=$PSPLINE_DIR \
    -DLIBSTELL_DIR=$LIBSTELL_DIR \
    -DSTRUMPACK_DIR=$STRUMPACK_DIR

# Build
make -j4
```

The executables `PB3D` and `POST` are created in the build directory.

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
| `STRUMPACK_DIR` | ~/Code/STRUMPACK-Dense-1.1.1 | Path to STRUMPACK-Dense 1.1.1 installation |

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
    -DLIBSTELL_DIR=$LIBSTELL_DIR \
    -DSTRUMPACK_DIR=$STRUMPACK_DIR
make -j4
```

## 8. HPC Cluster with Module System

On HPC systems with module-based environments:

```bash
# Load required modules (example)
module load gcc openmpi petsc slepc hdf5 netcdf-fortran scalapack

# Configure and build
mkdir build && cd build
cmake .. \
    -DPSPLINE_DIR=/path/to/pspline \
    -DLIBSTELL_DIR=/path/to/libstell \
    -DSTRUMPACK_DIR=/path/to/STRUMPACK-Dense-1.1.1
make -j4
```

## 9. Legacy Makefile

The original Makefile is preserved as `Makefile.legacy` for reference. To use it:

1. Copy it back: `cp Makefile.legacy Makefile`
2. Edit the paths at the top for your system
3. Run `make all`

The CMake-based build is recommended for new installations.

## 10. Sharing the setup

You can freeze the environment spec for collaborators:

```bash
spack env activate pb3d-env
spack concretize --fresh
spack env export > spack-lock.yaml
```

The `spack-lock.yaml` pins exact versions and hashes so other users reproduce the same build precisely.

## 11. Useful Spack commands

- `spack find` – list installed packages.
- `spack find --deps <spec>` – show dependency tree.
- `spack uninstall <spec>` – remove a package.
- `spack clean -a` – clear build caches if disk space is tight.
- `spack env location --view` – get path to environment view

For more background, see the official docs: <https://spack.readthedocs.io/>.
