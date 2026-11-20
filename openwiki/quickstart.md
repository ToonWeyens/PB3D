---
type: "Reference"
title: "PB3D — OpenWiki Quickstart"
description: "Entry point to the PB3D wiki: what PB3D and POST are, the repository layout, which page to read for which task (build, modules, physics, vacuum, testing, input, run sequence), and the facts to keep in mind before changing the code."
tags: [overview, quickstart, navigation, pb3d]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-d44494ef3e497fea81240ef8
    resource: repo://CMakeLists.txt
  - id: openwiki-source-47d02fca3524898d5aae2b3b
    resource: repo://LICENSE
  - id: openwiki-source-99e71edd5b71e76c5f03dd84
    resource: repo://Modules/input_ops.f90
  - id: openwiki-source-0e5296407709fa9254158dd8
    resource: repo://Modules/num_vars.f90
  - id: openwiki-source-5ebda9f0933cb8606c609e8f
    resource: repo://PB3D.f90
  - id: openwiki-source-8257526ded3898f9506f023a
    resource: repo://POST.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# PB3D — OpenWiki Quickstart

## What this repository is

<!-- openwiki: broken internal link [/README.md] link "/README.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
**PB3D** (Peeling Ballooning in 3D) is a Fortran HPC code that computes the linear ideal-MHD stability of peeling-ballooning modes in 3-D plasma equilibria, tokamaks and stellarators. It was written by Toon Weyens (2012 to now); the [`/README.md`](/README.md) holds the authors' affiliations and a changelog per version, which is the project's release notes.

PB3D reads an equilibrium computed elsewhere, from **VMEC** (a NetCDF `wout` file, through the LIBSTELL library) or **HELENA** (a text mapping file), sets up field-aligned perturbations, computes the vacuum response outside the plasma for a free boundary, and solves a generalised eigenvalue problem for ω² with **SLEPc/PETSc**; a negative ω² is an unstable mode. A second program, **POST**, reads PB3D's HDF5 output and writes derived quantities and plots as HDF5/XDMF for tools such as ParaView.

There is no server, API or database. The two programs `PB3D` and `POST` run from the command line on a namelist input file and an equilibrium file or PB3D output, usually under MPI.

## Repository layout

| Path | Contents |
|---|---|
| `PB3D.f90`, `POST.f90` | the two main programs |
| `Modules/` | the 50 Fortran modules with all physics, numerics, input/output and drivers |
| `Libraries/` | vendored FFT (`dfft.f`) and terminal output (`foul.f90`), and the build scripts and patch for PSPLINE and the minimal LIBSTELL |
| `include/` | preprocessor headers: the error macro `CHCKERR` (`PB3D_macros.h`) and the complex-number wrappers (`wrappers.h`) |
| `tests/` | the CTest suite: `unit/`, `fullstack/`, `regression/` (decks, scripts, checks), `fixtures/` (compressed equilibria), `reference/` (reference data and recipes), `vendor/` (test-drive) |
| `Examples/` | sample input decks and perturbation files |
| `Documentation/` | guides (`spack-setup.md`, `snellius-setup.md`, `testing.md`, `quality-and-ci.md`), the vacuum work (`vacuum_work_overview.md`, `vacuum_work_report.md`, `plan_vacuum_bugfix.md`, `w7x_check_vac_diag.md`), benchmark literature (`benchmarking/`) and figures; the Doxygen API documentation is at <https://pb3d.github.io/Doxygen/html/index.html> and built from `Doxyfile` |
| `scripts/` | Snellius build and job scripts with comparison helpers (`snellius/`), and the TERPSICHORE converter (`terpsichore/`) |
| `cmake/` | the `Find` modules for PSPLINE and LIBSTELL |
| `.github/` | the CI workflow, warning baselines and the apt mirror script |
| `CMakeLists.txt`, `spack.yaml`, `pyproject.toml`, `fortitude.toml`, `.fortlsrc` | build, dependency, Python, lint and language-server configuration |
| `Test/` | old standalone experiments with their own makefiles, not part of the build |
| `Makefile.legacy`, `ObjectList`, `PB3D.dep` | the old makefile build, not used by CMake |

## Which page to read

- **[Build & Dependencies](architecture/build-and-dependencies.md):** compiling with CMake, the build options, the libraries (PETSc/SLEPc 3.25 with complex scalars, HDF5, NetCDF, ScaLAPACK, optional STRUMPACK-Dense, PSPLINE, LIBSTELL) and how to build the ones Spack does not carry.
- **[Module Map](architecture/module-map.md):** what each module in `Modules/` does, in dependency order.
- **[Physics Domains](domain/physics-domains.md):** the equilibrium and its normalisation, perturbation modes, the radial discretisation and boundary conditions, the eigenvalue problem, Richardson extrapolation and POST.
- **[Vacuum Response](domain/vacuum-response.md):** the free-boundary boundary element method in its two styles, the diagonals, the response matrix, `check_vac`, the vacuum potential plot, its tests and its limits.
- **[Testing & Extension](testing-and-extension.md):** the test layers and how to run them, CI, exit status, versioning, and where to start a change.
- **[Input & Configuration](workflows/input-and-configuration.md):** the namelists, defaults and input checks, and how to add an option.
- **[Simulation Pipeline](workflows/simulation-pipeline.md):** what happens when PB3D and POST run, from start-up through the Richardson loop to the exit status.

## Facts to keep in mind

- **Two programs, one module library.** `PB3D` and `POST` link the same object library `pb3d_modules`, built from `Modules/*.f90`.
- **MPI from the first line.** Both programs call `start_MPI()` first and write parallel HDF5. Runs on more than one process need a distributed LU in PETSc, such as MUMPS.
- **The equilibrium decides much.** PB3D recognises VMEC or HELENA from the file (`eq_style` 1 or 2). VMEC runs use several field lines (`alpha_style = 2`) and the field-line 3-D vacuum; HELENA runs use one field line and the axisymmetric vacuum.
- **Boundary conditions.** `BC_style = 1 2` (zero at the inner end, free at the edge) is the default; the edge can also be fixed (1). No other values are accepted.
- **Tests exist; run them.** `ctest` in a full build runs the unit, full-stack and regression layers; the regression decks compare eigenvalues with recorded anchors. Extend the tests when you change the vacuum, the mode set-up or the matrix assembly.
- **Versions.** `prog_version` and `min_PB3D_version` in `num_vars.f90` (both 2.50) decide which PB3D output POST and a restart accept; change `min_PB3D_version` when the output layout changes, and add a changelog section to the README.
- **Companion tools.** Cluster run scripts and parameter scans for general use live in the separate `PB3D_tools` repository; `scripts/` here holds only the Snellius and TERPSICHORE scripts used during development.
<!-- openwiki: broken internal link [/LICENSE] link "/LICENSE" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- **Licence:** a BSD-style licence of the ITER Organization (2014-2019) with two extra conditions: the program keeps the name PB3D unless the original author agrees otherwise, and changes that lead to published or public work are to be sent back as a pull request (see [`/LICENSE`](/LICENSE)).
