!------------------------------------------------------------------------------!
!   Numerical utilities related to SLEPC (and PETSC) operations                !
!------------------------------------------------------------------------------!
! References:                                                                  !
! [1]   http://slepc.upv.es/documentation/current/docs/manualpages/EPS/        !
!       EPSComputeRelativeError.html                                           !
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
    public insert_block_mat, get_ghost_vec
#if ldebug
    public debug_insert_block_mat
#endif
    
    ! global variables
#if ldebug
    logical :: debug_insert_block_mat = .false.                                 ! plot debug information for insert_block_mat
#endif
    
contains
    ! Insert  a block pertaining  to 1 normal  point into a  matrix. Optionally,
    ! also set the Hermitian transpose and / or overwrite instead of add value.
    ! This  routine takes  into account  that if  the mode  numbers change  as a
    ! function of  the normal coordinate,  the local  blocks, which were  set up
    ! without  taking this  into  account,  have to  be  shifted. The  therefore
    ! missing information is omitted, and  the corresponding matrix elements are
    ! set to zero.
    ! This is  justifiable if there  are many more  mode couplings that  are not
    ! omitted  (i.e. when  n_mod_X is  large),  so that  it is  a small  effect,
    ! and/or  when  the coupling  between  these  modes is  already  negligible.
    ! Unfortunately, it  can be argued that  this is not generally  the case for
    ! the fast version (X_type = 2) since PV^1 and PV^2 are both ~ (nq-m), which
    ! indicates  that the  smallest elements  in the  corresponding off-diagonal
    ! blocks lie in the central columns for  PV^1 and in the central columns and
    ! rows for PV^2,  while the problematic elements lie farthest  away from the
    ! central columns  for PV^1 and  the central columns  and rows for  PV^2, so
    ! these are generally not small.
    ! The  coefficients  due to  normal  discretization  get increasingly  small
    ! for  higher orders,  though, which  compensates this  somewhat for  higher
    ! orders. Finally, an  alternative would be to set  the problematic elements
    ! not  to zero  but  equal  to the  closest  element,  but then  consistency
    ! would  be sacrified,  for  example  in the  energy  reconstruction of  the
    ! POST-processing: In  the energy  reconstruction the  exact same  terms are
    ! neglected as well.
    ! To conclude, for  fast X type 2,  the number of modes (n_mod_X)  has to be
    ! chosen high enough  compared with the variation of the  mode numbers; i.e.
    ! the variation of the safety factor or rotational transform.
    integer function insert_block_mat(block,mat,r_id,ind,n_sol,transp,&
        &overwrite) result(ierr)
        use X_vars, only: n_X, m_X, n_mod_X
        use num_vars, only: use_pol_flux_F
#if ldebug
        use input_utilities, only: pause_prog
#endif
        
        character(*), parameter :: rout_name = 'insert_block_mat'
        
        ! input / output
        PetscScalar :: block(:,:)                                               ! (n_mod x n_mod) block matrix for 1 normal point
        Mat, intent(inout) :: mat                                               ! matrix in which to insert block
        PetscInt, intent(in) :: r_id                                            ! normal position of corresponding V^0 (starting at 0)
        PetscInt, intent(in) :: ind(2)                                          ! 2D index in matrix, relative to r_id
        PetscInt, intent(in) :: n_sol                                           ! number of grid points of solution grid
        PetscBool, intent(in), optional :: transp                               ! also set Hermitian transpose
        PetscBool, intent(in), optional :: overwrite                            ! overwrite
        
        ! local variables
        PetscInt, allocatable :: loc_k(:), loc_m(:)                             ! the locations at which to add the blocks to the matrices
        PetscInt :: kd                                                          ! counter
        PetscInt :: k, m                                                        ! counters
        PetscInt :: k_loc, m_loc                                                ! local k and m
        PetscBool :: transp_loc                                                 ! local copy of transp
        PetscBool :: overwrite_loc                                              ! local copy of overwrite
        character(len=max_str_ln) :: err_msg                                    ! error message
        PetscInt :: operation                                                   ! either ADD_VALUES or INSERT_VALUES
        PetscScalar, allocatable :: block_loc(:,:)                              ! local block, possibly shifted from block
        
        ! initialize ierr
        ierr = 0
        
        ! set up local transp and overwrite
        transp_loc = .false.
        if (present(transp)) transp_loc = transp
        overwrite_loc = .false.
        if (present(overwrite)) overwrite_loc = overwrite
        
        ! set up local k and m
        allocate(loc_k(size(block,1)))
        allocate(loc_m(size(block,2)))
        loc_k = [(kd, kd = 0,n_mod_X-1)] + (r_id+ind(1))*n_mod_X
        loc_m = [(kd, kd = 0,n_mod_X-1)] + (r_id+ind(2))*n_mod_X
        
        ! set operation
        if (overwrite_loc) then
            operation = INSERT_VALUES
        else
            operation = ADD_VALUES
        end if
        
#if ldebug
        ! user output
        if (debug_insert_block_mat) then
            call writo('>>> at (k,m) = '//trim(i2str(r_id))//' + ('//&
                &trim(i2str(ind(1)))//','//trim(i2str(ind(2)))//'):')
            call lvl_ud(1)
        end if
#endif
        
        ! only set values if within matrix range
        if (minval(r_id+ind).ge.0 .and. maxval(r_id+ind).lt.n_sol) then
            ! set error message
            err_msg = 'Couldn''t add values to matrix'
            
            ! set up local block, possibly translated for fast PB3D
            allocate(block_loc(size(block,1),size(block,2)))
            block_loc = 0._dp
            do m = 1,n_mod_X
                do k = 1,n_mod_X
                    if (use_pol_flux_F) then
                        k_loc = k + m_X(r_id+1,k) - m_X(r_id+1+ind(1),k)
                        m_loc = m + m_X(r_id+1,k) - m_X(r_id+1+ind(2),k)
                    else
                        k_loc = k + n_X(r_id+1,k) - n_X(r_id+1+ind(1),k)
                        m_loc = m + n_X(r_id+1,k) - n_X(r_id+1+ind(2),k)
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
                call writo('following local block is going to be added: ')
                call writo('Re =')
                call print_ar_2(rp(block_loc))
                call writo('Im =')
                call print_ar_2(ip(block_loc))
                if (transp_loc) call writo('as well as at the transposed place')
                if (overwrite_loc) then
                    call writo('(with operation INSERT_VALUES)')
                else
                    call writo('(with operation INSERT_VALUES)')
                end if
            endif
#endif
            
            ! set values
            call MatSetValues(mat,n_mod_X,loc_k,n_mod_X,loc_m,&
                &transpose(block_loc),operation,ierr)
            CHCKERR(err_msg)
            
            ! set transpose if wanted
            if (transp_loc) then
                call MatSetValues(mat,n_mod_X,loc_m,n_mod_X,loc_k,&
                    &conjg(block_loc),operation,ierr)
                CHCKERR(err_msg)
            end if
            
            ! clean up
            deallocate(block_loc)
        else
#if ldebug
            if (debug_insert_block_mat) &
                &call writo('Due to out of range no local block is added')
#endif
        end if
        
        ! clean up
        deallocate(loc_k,loc_m)
        
#if ldebug
        if (debug_insert_block_mat) then
            call pause_prog()
            call lvl_ud(-1)
        endif
#endif
    end function insert_block_mat
    
    ! gets ghost region of a Petsc Vector.
    !!!!! The ghost regions at the end of the total range are set to zero.
    ! Note: Every process should have at  least a number of normal poins greater
    ! than n_ghost. For performance, this is not checked.
    integer function get_ghost_vec(V_ptr,V_ghost,n_ghost) result(ierr)
        use num_vars, only: rank, n_procs
        use MPI
        
        character(*), parameter :: rout_name = 'get_ghost_vec'
        
        ! input / output
        PetscScalar, intent(in), pointer :: V_ptr(:)                            ! pointer to vector of which ghost zones are to be retrieved
        PetscScalar, intent(inout), allocatable :: V_ghost(:,:)                 ! ghost zone for previous (1) and next process (2)
        integer, intent(in) :: n_ghost                                          ! size of ghost region
        
        ! local variables
        integer :: n_loc                                                        ! local size
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! ghost regions only make sense if there is more than 1 process
        if (n_procs.gt.1) then
            ! set up ghost regions
            n_loc = size(V_ptr)
            allocate(V_ghost(n_ghost,2))
            
            ! communicate with previous proc
            if (rank.eq.0) then                                                 ! first rank only receives
                call MPI_Recv(V_ghost(:,2),n_ghost,MPI_DOUBLE_COMPLEX,&
                    &rank+1,rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else if (rank+1.eq.n_procs) then                                    ! last rank only sends
                !V_ghost(:,2) = 0._dp
                V_ghost(:,2) = huge(1._dp)
                call MPI_Send(V_ptr(1:n_ghost),n_ghost,MPI_DOUBLE_COMPLEX,&
                    &rank-1,rank,MPI_Comm_world,ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(V_ptr(1:n_ghost),n_ghost,&
                    &MPI_DOUBLE_COMPLEX,rank-1,rank,&
                    &V_ghost(:,2),n_ghost,MPI_DOUBLE_COMPLEX,rank+1,rank+1,&
                    &MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
            
            ! communicate with next proc
            if (rank.eq.0) then                                                 ! first rank only sends
                !V_ghost(:,1) = 0._dp
                V_ghost(:,1) = huge(1._dp)
                call MPI_Send(V_ptr(n_loc-n_ghost:n_loc-1),n_ghost,&
                    &MPI_DOUBLE_COMPLEX,rank+1,rank,MPI_Comm_world,ierr)
                CHCKERR('Failed to send')
            else if (rank+1.eq.n_procs) then                                    ! last rank only receives
                call MPI_Recv(V_ghost(:,1),n_ghost,MPI_DOUBLE_COMPLEX,&
                    &rank-1,rank-1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(V_ptr(n_loc-n_ghost:n_loc-1),n_ghost,&
                    &MPI_DOUBLE_COMPLEX,rank+1,rank,&
                    &V_ghost(:,1),n_ghost,MPI_DOUBLE_COMPLEX,rank-1,rank-1,&
                    &MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_vec
end module SLEPC_utilities
