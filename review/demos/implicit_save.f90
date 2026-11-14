! Area 0a: the PETSc port fixed an implicit SAVE in SLEPC_ops (remove_next_loc).
!
! At 2.47 the local variable was declared as
!     PetscBool :: remove_next_loc = .false.
! An initialisation in a declaration gives the variable the SAVE attribute: it
! is initialised once, at program start, and keeps its value between calls.
! After one call with remove_next = .true., every later call without the
! argument still sees .true.. The port declares it without initialisation and
! sets it to .false. at the start of every call.
!
! Build and run:  gfortran implicit_save.f90 -o implicit_save && ./implicit_save
program implicit_save
    implicit none
    call as_247(.true.)
    call as_247()
    call as_247()
    print *
    call as_now(.true.)
    call as_now()
    call as_now()
contains
    subroutine as_247(remove_next)
        logical, intent(in), optional :: remove_next
        logical :: remove_next_loc = .false.        ! implies SAVE
        if (present(remove_next)) remove_next_loc = remove_next
        print '(a,l1,a,l1)', '2.47: argument present = ', present(remove_next), &
            &', remove_next_loc = ', remove_next_loc
    end subroutine as_247
    subroutine as_now(remove_next)
        logical, intent(in), optional :: remove_next
        logical :: remove_next_loc
        remove_next_loc = .false.
        if (present(remove_next)) remove_next_loc = remove_next
        print '(a,l1,a,l1)', 'now:  argument present = ', present(remove_next), &
            &', remove_next_loc = ', remove_next_loc
    end subroutine as_now
end program implicit_save
