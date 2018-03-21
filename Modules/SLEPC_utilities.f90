!------------------------------------------------------------------------------!
!> Numerical utilities related to SLEPC (and PETSC) operations.
!------------------------------------------------------------------------------!
!> \see
!! References:
!! \cite Hernandez2005slepc
!------------------------------------------------------------------------------!
module SLEPC_utilities
#include <PB3D_macros.h>
#include <wrappers.h>
! for slepc 3.6.0:
#include <slepc/finclude/slepcepsdef.h>
! for slepc 3.5.3:
!#include <finclude/slepcepsdef.h>
    use str_utilities
    use output_ops
    use messages
    use slepceps
    use num_vars, only: iu, dp, max_str_ln

    implicit none
    private
    public insert_block_mat
#if ldebug
    public debug_insert_block_mat
#endif
    
    ! global variables
#if ldebug
    logical :: debug_insert_block_mat = .false.                                 !< plot debug information for insert_block_mat() \ldebug
#endif
    
contains
    !> Insert a block pertaining to 1 normal point into a matrix.
    !!
    !! Optionally, also set  the Hermitian transpose and /  or overwrite instead
    !! of add value.
    !!
    !! This routine  takes into  account that  if the mode  numbers change  as a
    !! function of  the normal coordinate, the  local blocks, which were  set up
    !! without  taking this  into account,  have  to be  shifted. The  therefore
    !! missing information is omitted, and the corresponding matrix elements are
    !! set to zero.
    !!
    !! This is  justifiable if there are  many more mode couplings  that are not
    !! omitted (i.e. when  \c n_mod_X is large),  so that it is  a small effect,
    !! and/or when the coupling between these modes is already negligible.
    !!
    !! Unfortunately, it can  be argued that this is not  generally the case for
    !! the  fast  version (\c  X_style  =  2) since  \f$\widetilde{PV}^1\f$  and
    !! \f$\widetilde{PV}^2\f$ are  both \f$\sim (nq-m)\f$, which  indicates that
    !! the smallest elements in the corresponding off-diagonal blocks lie in the
    !! central columns for \f$\widetilde{PV}^1\f$ and in the central columns and
    !! rows  for  \f$\widetilde{PV}^2\f$,  while the  problematic  elements  lie
    !! farthest away from the central columns for \f$\widetilde{PV}^1\f$ and the
    !! central  columns  and  rows  for  \f$\widetilde{PV}^2\f$,  so  these  are
    !! generally not small.
    !!
    !! The coefficients due to normal  discretization get increasingly small for
    !! higher orders, though, which compensates this for higher orders.
    !!
    !! Finally,  an  alternative  would  be  to  set  the  problematic  elements
    !! not  to zero  but  equal to  the closest  element,  but then  consistency
    !! would  be sacrified,  for example  in  the energy  reconstruction of  the
    !! POST-processing: In  the energy reconstruction  the exact same  terms are
    !! neglected as well.
    !!
    !! To conclude, for fast \c X_style 2,  the number of modes (\c n_mod_X) has
    !! to be chosen high enough compared with the variation of the mode numbers;
    !! i.e. the variation of the safety factor or rotational transform.
    !!
    !! \return ierr
    integer function insert_block_mat(mds,block,mat,r_id,ind,n_r,transp,&
        &overwrite,ind_insert) result(ierr)
        use X_vars, only: modes_type, &
            &n_mod_X
        use num_vars, only: use_pol_flux_F
#if ldebug
        use num_vars, only: rank
        use input_utilities, only: pause_prog
#endif
        
        character(*), parameter :: rout_name = 'insert_block_mat'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        PetscScalar, intent(in) :: block(:,:)                                   !< (\cn_mod x \cn_mod) block matrix for 1 normal point
        Mat, intent(inout) :: mat                                               !< matrix in which to insert block
        PetscInt, intent(in) :: r_id                                            !< normal position of corresponding \f$\widetilde{V}^0\f$ (starting at 0)
        PetscInt, intent(in) :: ind(2)                                          !< 2D index in matrix, relative to \cr_id
        PetscInt, intent(in) :: n_r                                             !< number of grid points of solution grid
        PetscBool, intent(in), optional :: transp                               !< also set Hermitian transpose
        PetscBool, intent(in), optional :: overwrite                            !< overwrite
        PetscBool, intent(in), optional :: ind_insert                           !< individual insert, only important for debugging
        
        ! local variables
        PetscInt :: k, m                                                        ! counters
        PetscInt :: k_loc, m_loc                                                ! local k and m
        PetscBool :: transp_loc                                                 ! local copy of transp
        PetscBool :: overwrite_loc                                              ! local copy of overwrite
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: operation                                                   ! either ADD_VALUES or INSERT_VALUES
        PetscScalar, allocatable :: block_loc(:,:)                              ! local block, possibly shifted from block
        PetscBool :: ind_insert_loc                                             ! local ind_insert
        
        ! initialize ierr
        ierr = 0
        
        ! set up local transp, overwrite and ind_insert
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        overwrite_loc = .false.
        if (present(overwrite)) overwrite_loc = overwrite
        ind_insert_loc = .false.
        if (present(ind_insert)) ind_insert_loc = ind_insert
        
        ! set operation
        if (overwrite_loc) then
            operation = INSERT_VALUES
        else
            operation = ADD_VALUES
        end if
        
#if ldebug
        ! user output
        if (debug_insert_block_mat) then
            if (.not.ind_insert_loc) call sleep(rank)
            call writo('>>> rank '//trim(i2str(rank))//&
                &': at (k,m) = '//trim(i2str(r_id))//' + ('//&
                &trim(i2str(ind(1)))//','//trim(i2str(ind(2)))//'):',&
                &persistent=.true.)
            call lvl_ud(1)
        end if
#endif
        
        ! only set values if within matrix range
        if (minval(r_id+ind).ge.0 .and. maxval(r_id+ind).lt.n_r) then
            ! set error message
            err_msg = 'Couldn''t add values to matrix'
            
            ! set up local block, possibly translated for fast PB3D
            ! Note: For BC_style 3, there is  a grid extension beyond n_r, which
            ! is why there is a 'min(..,size(n_X,1))' here.
            allocate(block_loc(size(block,1),size(block,2)))
            block_loc = 0._dp
            do m = 1,n_mod_X
                do k = 1,n_mod_X
                    if (use_pol_flux_F) then
                        k_loc = k + mds%m(r_id+1,k) - &
                            &mds%m(min(r_id+1+ind(1),size(mds%m,1)),k)
                        m_loc = m + mds%m(r_id+1,k) - &
                            &mds%m(min(r_id+1+ind(2),size(mds%m,1)),k)
                    else
                        k_loc = k + mds%n(r_id+1,k) - &
                            &mds%n(min(r_id+1+ind(1),size(mds%n,1)),k)
                        m_loc = m + mds%n(r_id+1,k) - &
                            &mds%n(min(r_id+1+ind(2),size(mds%n,1)),k)
                    end if
                    if (k_loc.ge.1 .and. m_loc.ge.1 .and. &
                        &k_loc.le.n_mod_X .and. m_loc.le.n_mod_X) &
                        block_loc(k_loc,m_loc) = block(k,m)
                end do
            end do
            
            !write(*,*) 'block'
            !do k = 1,n_mod_X
                !write(*,*) k, block(k,:)
            !end do
            !write(*,*) 'block_loc'
            !do k = 1,n_mod_X
                !write(*,*) k, block_loc(k,:)
            !end do
            
#if ldebug
            if (debug_insert_block_mat) then
                call writo('following local block is going to be added by &
                    &rank '//trim(i2str(rank))//':',persistent=.true.)
                call writo('Re =',persistent=.true.)
                call print_ar_2(rp(block_loc))
                call writo('Im =',persistent=.true.)
                call print_ar_2(ip(block_loc))
                if (transp_loc) call writo('as well as at the transposed place')
                if (overwrite_loc) then
                    call writo('(with operation INSERT_VALUES)',&
                        &persistent=.true.)
                else
                    call writo('(with operation ADD_VALUES)',&
                        &persistent=.true.)
                end if
            endif
#endif
            
            ! set values
            call MatSetValuesBlocked(mat,1,r_id+ind(1),1,r_id+ind(2),&
                &transpose(block_loc),operation,ierr)
            CHCKERR(err_msg)
            
            ! set transpose if wanted
            if (transp_loc) then
                call MatSetValuesBlocked(mat,1,r_id+ind(2),1,r_id+ind(1),&
                    &conjg(block_loc),operation,ierr)
                CHCKERR(err_msg)
            end if
            
            ! clean up
            deallocate(block_loc)
        else
#if ldebug
            if (debug_insert_block_mat) &
                &call writo('Due to out of range no local block is added',&
                &persistent=.true.)
#endif
        end if
        
#if ldebug
        if (debug_insert_block_mat) then
            call pause_prog(ind_insert_loc)
            call lvl_ud(-1)
        endif
#endif
    end function insert_block_mat
end module SLEPC_utilities
