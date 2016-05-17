- see http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html for other ways to solve
  http://slepc.upv.es/documentation/slepc.pdf#subsection.3.4.2 explains this.
- See http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex14f.F.html for basic example of using matcreateshell.
- See http://www.mcs.anl.gov/petsc/petsc-3.6/src/snes/examples/tutorials/ex5f90.F.html for instructions on how to use context in FORTRAN.
- See http://slepc.upv.es/documentation/slepc.pdf#section.8.2 in SLEPC manual to read about shell matrices.
- See http://slepc.upv.es/documentation/slepc.pdf#section.3.4 in SLEPC manual to read about spectral transforms.
- lu preconditioning is difficult for shell matrices
- Jacobi preconditioning only needs the Diagonal, so it is easy to use, but iterative.
