# Vendored: test-drive

- Project: <https://github.com/fortran-lang/test-drive>
- Version: v0.5.0 (file `src/testdrive.F90`, unmodified)
- License: Apache-2.0 OR MIT (see `LICENSE-MIT`)

test-drive is the fortran-lang community's light-weight unit testing
framework (also used by fortran stdlib). It is vendored here as a single
source file so that PB3D's test suite builds offline and on HPC systems
without package-manager access.

To upgrade: replace `testdrive.F90` with a newer tagged version and update
this note.
