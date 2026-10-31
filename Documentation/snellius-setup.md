# PB3D on Snellius

How PB3D is built and run on the Dutch national supercomputer Snellius
(SURF, `snellius.surf.nl`), with the scripts under `scripts/snellius/`.
Set up 2026-10-27/28 on login node `int5`, Rome compute nodes, partition
`rome`. In the commands below `snellius` stands for
`tweyens@snellius.surf.nl`.

Everything lives under `~/pb3d` on Snellius. No module, library or file
outside that directory was touched, and `rm -rf ~/pb3d` removes all of it.

## Layout on Snellius

| Path | Content |
|---|---|
| `~/pb3d/env.sh` | modules and environment (source it in every shell and job) |
| `~/pb3d/scripts/` | the scripts of this directory (`scripts/snellius/`) |
| `~/pb3d/src/` | clean source tree (`git archive HEAD` of the worktree) |
| `~/pb3d/deps/petsc/` | PETSc 3.25.3, complex scalars, with MUMPS/ParMETIS/METIS, prefix install |
| `~/pb3d/deps/slepc/` | SLEPc 3.25.1, prefix install |
| `~/pb3d/deps/pspline/` | PSPLINE checkout (`install/` holds `lib/libpspline.a` + `.mod`) |
| `~/pb3d/deps/libstell-min/` | minimal LIBSTELL (`read_wout_mod` and dependencies) |
| `~/pb3d/deps/stellopt/LIBSTELL/Sources/` | the STELLOPT LIBSTELL sources (rsynced, 8.5 MB) |
| `~/pb3d/deps/venv/` | Python venv holding CMake 4.4.4 (see "Known points") |
| `~/pb3d/fixtures/cbm18a_extended_vac.12` | the HELENA equilibrium fixture (25 MB) |
| `~/pb3d/build/`, `~/pb3d/build-debug/` | Release and Debug CMake builds |
| `~/pb3d/runs/<job>-<jobid>/np<N>/` | one directory per run |
| `~/pb3d/logs/` | build logs and Slurm output (`%x-%j.out`) |

The PETSc and SLEPc build trees and tarballs are deleted after `make
install`; only the install prefixes stay (all of `~/pb3d/deps` is 241 MB, of which the CMake venv 49 MB).

## Modules

`env.sh` does `module purge` and loads, from the 2025 software stack:

- `foss/2025a`: GCC 14.2.0, OpenMPI 5.0.7, FlexiBLAS 3.4.5 over OpenBLAS
  0.3.29, ScaLAPACK 2.2.2
- `HDF5/1.14.6-gompi-2025a` (parallel, with Fortran bindings)
- `netCDF-Fortran/4.6.2-gompi-2025a`
- `CMake/3.31.3-GCCcore-14.2.0`
- `pkgconf/2.3.0-GCCcore-14.2.0`

and then puts `~/pb3d/deps/venv/bin` (CMake 4.4.4) in front of the PATH.

The site `PETSc/3.23.5-foss-2025a` module is built with real scalars and
cannot be used: PB3D needs complex PETSc, so it is built from source. The
site `MUMPS/5.8.1-foss-2025a-metis` module is not used either; PETSc
downloads and builds its own MUMPS (simpler to configure, tested
combination).

## Build recipe

```bash
# on the Mac: source tree, fixture, STELLOPT sources, scripts
cd <PB3D worktree>
ssh snellius 'mkdir -p ~/pb3d/{src,fixtures,scripts,deps/stellopt/LIBSTELL/Sources}'
git archive HEAD | ssh snellius 'tar -x -C ~/pb3d/src'
scp Libraries/pspline-r8genxpkg-nx3.patch snellius:~/pb3d/src/Libraries/   # if not yet committed
rsync -az ~/Code/PB3D-fixtures/cbm18a_extended_vac.12 snellius:~/pb3d/fixtures/
rsync -az ~/Code/stellopt/LIBSTELL/Sources/ snellius:~/pb3d/deps/stellopt/LIBSTELL/Sources/
rsync -az scripts/snellius/ snellius:~/pb3d/scripts/
scp scripts/snellius/env.sh snellius:~/pb3d/env.sh

# on Snellius (login node; nice, at most 16 jobs)
~/pb3d/scripts/build_deps.sh all        # petsc | slepc | pspline | libstell | all
~/pb3d/scripts/build_pb3d.sh all        # release | debug | all
# or everything above plus the jobs of the next section, unattended:
nohup ~/pb3d/scripts/pipeline.sh > ~/pb3d/logs/pipeline.log 2>&1 &
```

`build_deps.sh` downloads anonymously (generic User-Agent, no identity in
any request): the PETSc release tarball from `web.cels.anl.gov`, the SLEPc
tarball from `slepc.upv.es`, and `git clone` of
`github.com/ilonster/pspline`. STELLOPT is not cloned; only
`LIBSTELL/Sources` is rsynced from the Mac.

PETSc configure line (in `build_deps.sh`):
`--with-scalar-type=complex --with-precision=double --with-debugging=0
--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90
--with-blaslapack-lib="-L$EBROOTFLEXIBLAS/lib -lflexiblas"
--with-scalapack-lib="-L$EBROOTSCALAPACK/lib -lscalapack"
--download-mumps --download-parmetis --download-metis
--with-fortran-bindings=1 --with-shared-libraries=1 --with-x=0
COPTFLAGS=-O2 CXXOPTFLAGS=-O2 FOPTFLAGS=-O2 --prefix=~/pb3d/deps/petsc`.
SLEPc: `./configure --prefix=~/pb3d/deps/slepc` against that PETSc. Both
leave `lib/pkgconfig/{PETSc,SLEPc}.pc`, which CMake finds through
`PKG_CONFIG_PATH` from `env.sh`.

PSPLINE: `build_pspline_linux.sh` is `Libraries/build_pspline.sh` without
the Homebrew `PATH` line, taking the checkout as its first argument; it
applies `Libraries/pspline-r8genxpkg-nx3.patch` when `build_deps.sh` finds
it in `~/pb3d/src/Libraries/`. LIBSTELL: `build_libstell_min_linux.sh` is
`Libraries/build_libstell_min.sh` with the Homebrew line removed and the
defaults pointing at `$EBROOTNETCDFMINFORTRAN` and `~/pb3d/deps/stellopt`.

CMake configure (in `build_pb3d.sh`):
`cmake -S ~/pb3d/src -B ~/pb3d/build -DCMAKE_BUILD_TYPE=Release
-DPSPLINE_DIR=~/pb3d/deps/pspline/install
-DLIBSTELL_DIR=~/pb3d/deps/libstell-min -DPB3D_FIXTURE_DIR=~/pb3d/fixtures`;
the Debug variant adds `-DCMAKE_BUILD_TYPE=Debug -DPB3D_ENABLE_DEBUG=ON`
in `~/pb3d/build-debug`. STRUMPACK is not built; the vacuum uses the
ScaLAPACK `pdgesv` path.

Build times on the login node (`nice -n 10`, `-j16`): PETSc 21 min
without MUMPS, 35 min with MUMPS 5.8.2, ParMETIS and METIS downloaded and built by its configure; SLEPc 4.7
min, PSPLINE 2 min, LIBSTELL under 1 min, PB3D Release 144 s, PB3D Debug
116 s. The CMake configure of PB3D takes 10 s.

## Running a job

Rules followed here (and to follow): one node per job, partition `rome`,
`--time` of 15 minutes or less (30 for `mem_run.sbatch`), output under
`~/pb3d/logs/%x-%j.out` (the scripts write `/home/%u/...`, which Slurm
fills with the user name), nothing left queued or running afterwards.
Several jobs may run at the same time (Toon, 30 October 2026); the night
of 30 October had up to six in the queue.

```bash
source ~/pb3d/env.sh
D=~/pb3d/src/tests/regression
# run_deck.sbatch <deck> [np-list] [bc-style] [repeat] [exe] [keep] [eq-file]
sbatch ~/pb3d/scripts/run_deck.sbatch $D/input_cbm18a_fixed_bnd "0 1 2 4"        # 0 = singleton, else mpirun -np N
sbatch ~/pb3d/scripts/run_deck.sbatch $D/input_cbm18a_free_bnd "0 1 2 4" "1 2"   # BC_style override
sbatch ~/pb3d/scripts/run_deck.sbatch $D/input_cbm18a_fixed_bnd "0" "" 10        # 10 repeats: crash statistics
sbatch ~/pb3d/scripts/run_deck.sbatch ~/pb3d/decks/input_twin_fixed "0" "" 1 "" "" \
    ~/pb3d/fixtures/wout_circ_eps033_q104.nc                                      # another equilibrium (VMEC: keep the wout_*.nc name)
# mem_run.sbatch <exe> <deck> <eq-file>: one run under /usr/bin/time -v, prints wall time and peak memory
sbatch ~/pb3d/scripts/mem_run.sbatch ~/pb3d/build/PB3D ~/pb3d/decks/input_w7x_fixed ~/pb3d/fixtures/wout_w7x.nc
sbatch --job-name=pb3d-ctest --partition=rome --nodes=1 --ntasks=16 --time=00:15:00 \
    --output=$HOME/pb3d/logs/%x-%j.out \
    --wrap="source ~/pb3d/env.sh; cd ~/pb3d/build && ctest --output-on-failure --timeout 600"
```

`run_deck.sbatch` (16 tasks, 15 minutes) runs the deck once per process
count in `~/pb3d/runs/<job>-<jobid>/np<N>/`, with the deck copied to
`input` and the equilibrium symlinked as `eq` (PB3D limits file paths to
120 characters), and prints `EV np=<N> ... <eigenvalue>` lines (the first
accepted eigenvalue from `PB3D_out_EV_R_1.txt`). Each run writes a
`PB3D_out.h5` of about 270 MB; the script deletes it after reading the
eigenvalue unless the sixth argument is `keep` (the first day's runs had
piled up 19 GB before this was added). The options are
positional on purpose: Snellius submits with the environment of the
submitting shell dropped, so `BC_STYLE=... sbatch` or `--export=ALL,...`
does not reach the job (the first attempt here ran the "BC style 2" and
"10 repeats" jobs with the deck defaults). `pipeline.sh` runs the whole
sequence unattended: dependencies, the two builds, then the jobs one after
the other with `sbatch --wait`.

`mem_run.sbatch` runs one executable on one deck under `/usr/bin/time -v`
and prints one `MEM ...` line with the wall time, the peak resident memory
and the eigenvalue; it was used for the before/after table of step 9 of
issue #4. (`sacct`'s MaxRSS is empty or wrong for jobs shorter than about
a minute, because the accounting samples every 30 s.) The `compare_*.py`
helpers (h5py, numpy) read the `PB3D_out.h5` files kept with `keep`:
`compare_eq.py` compares the equilibrium and grid datasets of two outputs,
`compare_np.py` finds the radial indices where the integrated perturbation
matrices differ between process counts, `compare_petsc_mats.py` compares
PETSc binary matrix dumps, `compare_res.py` and `compare_res_series.py`
look at the stored 3-D vacuum response (diagonal ratios against a
reference, eigenvalues of the Hermitian part, convergence of a series of
runs, given as output files or run directories), and `checkvac_table.py`
turns the `CHECKVAC` lines of the run logs (`run.log` in the run
directory, written with `check_vac = .true.`) into one table row per run.
They were written for the investigations of 30 and 31 October 2026 and
are kept as the record of how those tables were made.

A run of the fixed-boundary deck takes 4-5 s on one Rome core, the
free-boundary deck 8-10 s; a whole `run_deck.sbatch` job with four process
counts finishes in under a minute, `ctest` of the Release build in 70 s.

## Results (2026-10-28)

Release build, PETSc 3.25.3 complex with MUMPS 5.8.2, SLEPc 3.25.1, patched
PSPLINE, one Rome node, partition `rome`. First accepted eigenvalue
\(\omega^2/\omega_A^2\) of `PB3D_out_EV_R_1.txt`; "0" is a singleton run
without `mpirun`.

| Deck | procs | eigenvalue | Mac (gfortran 16, PETSc 3.25.3) |
|---|---|---|---|
| `input_cbm18a_fixed_bnd` (BC_style 1 1) | 0 | -4.0735999821208513E-02 | -4.0735999821E-02 |
| | 1 | -4.0735999821208513E-02 | |
| | 2 | -4.0735232666401666E-02 | |
| | 4 | -4.0735497841722867E-02 | |
| `input_cbm18a_free_bnd` (BC_style 1 4) | 0 | -4.2348776924512666E-02 | -4.23487769E-02 |
| | 1 | -4.2348776924512666E-02 | |
| | 2 | -4.2348775506777903E-02 | |
| | 4 | -4.2344229173739534E-02 | |
| `input_cbm18a_free_bnd`, BC_style 1 2 | 0 | -4.2348778421032129E-02 | -4.23487784E-02 |
| | 1 | -4.2348778421032129E-02 | |
| | 2 | -4.2348778559244016E-02 | |
| | 4 | -4.2344214857731988E-02 | |

The singleton and 1-process values agree with the Mac to every printed
digit. With 2 and 4 processes the eigenvalues move by 1.9e-5 and 1.2e-5
(fixed), 3.3e-8 and 1.1e-4 (free, style 4), 3.3e-9 and 1.1e-4 (free,
style 2) relative; that is well inside the regression tolerance (10 %)
but far above the SLEPc tolerance (1e-6), so the distributed problem is
not the serial one to round-off. The 2-process fixed value came out
bit-identical in a second job, so the shift is systematic, not noise.
Found and fixed on 30 October 2026 (issue #4, commit 798fe50):
`redistribute_var` filled the overlap between two process ranges from the
lower rank's ghost region; since then the process counts agree to 1.5e-8,
the reproducibility of the eigenvalue itself (see
`Documentation/plan_vacuum_bugfix.md`, "Night of 30 October 2026").

Crash statistics, fixed deck, singleton, ten launches in a row: 10 of 10
produced the eigenvalue file with the unpatched PSPLINE and 10 of 10 with
the patched one, every value identical to the table. Runs with 2 and 4
processes before the MUMPS rebuild: 0 of 12 (see "Known points").

`ctest` (Release, 21 tests: 7 unit, 2 smoke, 10 full-stack including the
`np2`/`np4` ScaLAPACK rows and `read_HEL`, 2 regression): 21/21 passed in 69.5 s.
`ctest -L 'unit|smoke|fullstack'` on the Debug build (18 tests):
18/18 passed in 83.8 s. Before the MUMPS rebuild the same two runs gave 21/21
in 69.6 s and 18/18 in 84.0 s.

## Known points

- **CMake 3.31 loops on PETSc 3.25's Fortran headers.** With the site
  `CMake/3.31.3` (and the system 3.31.8), `cmake --build` never starts
  compiling: `cmake -E cmake_depends` for the `pb3d_modules` target runs
  for a quarter of an hour at 13 GB peak memory, opening
  `petsc/finclude/petscmath.h` and `petscsys.h` over and over. The two
  headers include each other under `#if !defined(...)` guards, which CMake
  3.31's own Fortran preprocessor does not terminate on. CMake 4.4.4 (the
  Mac builds use 4.3) scans the tree in a second, so `env.sh` puts the
  venv's CMake first. The Makefile generator is kept (same as on the Mac).
- **More than one MPI process needs MUMPS.** PB3D's eigensolve is
  shift-invert (`STSINVERT`), whose default KSP is a direct LU of the
  shifted operator. PETSc's built-in LU is sequential only, so with 2 or 4
  processes every run aborts in `EPSSolve` with "Could not locate a solver
  type for factorization type LU and matrix type mpibaij" (PB3D reports it
  as error code 92 from `get_solution`). The Mac's Spack PETSc is
  `~mumps ~superlu-dist`, and `mpirun` crashes on the Mac (Homebrew OpenMPI
  5), so this was never visible there. The recipe builds PETSc with
  `--download-mumps --download-parmetis --download-metis`; PETSc then picks
  MUMPS by itself for distributed LU, no PB3D change needed.
- **Slurm drops the submitting environment.** `BC_STYLE=... sbatch` and
  `sbatch --export=ALL,...` do not reach the job on Snellius; the job
  options of `run_deck.sbatch` are therefore positional.
- **PSPLINE `r8genxpkg` overrun** (`Libraries/pspline-r8genxpkg-nx3.patch`):
  invisible on Linux. Ten singleton runs of the fixed deck with the
  unpatched library all produced the eigenvalue file with the identical
  value; the patched library gives the same value (table above).
- `ld` warns once that `grid_utilities.f90.o` "requires executable stack"
  (its `.note.GNU-stack` section is marked executable; gfortran does this
  when an internal procedure is passed as an actual argument and needs a
  trampoline). Harmless for these runs; `-Wl,-z,noexecstack` would hide it.
- A singleton PB3D run (`PB3D input eq` without `mpirun`) works on the
  compute nodes; OpenMPI 5.0.7 `mpirun -np N` inside the Slurm allocation
  works too (no `srun` needed).
- ScaLAPACK solves the vacuum boundary-element system (`pdgesv`);
  STRUMPACK-Dense is not built.
