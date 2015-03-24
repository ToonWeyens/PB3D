!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi

    implicit none
    private
    public calc_real_X
    
    ! interfaces
    interface calc_real_X
        module procedure calc_real_X_arr, calc_real_X_ind
    end interface
    !interface calc_real_U
        !module procedure calc_real_U_arr, calc_real_U_ind
    !end interface

contains
    ! calculates the  normal component  of the  perturbation, or  optionally the
    ! parallel derivative. The input is given for a grid in Flux coordinates and
    ! for a range in normalized time (1 corresponds to one period).
    integer function calc_real_X_arr(eq,grid_X,X,X_id,time,X_F,deriv) &
        &result(ierr)                                                           ! (time) array version
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        use num_vars, only: use_pol_flux_X, grp_rank, grp_n_procs
        
        character(*), parameter :: rout_name = 'calc_real_X_arr'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: X_F(:,:,:,:)                                 ! normal component of perturbation
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        complex(dp), allocatable :: deriv_mult_factor(:)                        ! multiplicative factor due to parallel derivative
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(grid_X%theta_F,1)
        n_zeta = size(grid_X%theta_F,2)
        n_r = size(grid_X%theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(grid_X%zeta_F,1) .or. &
            &n_zeta.ne.size(grid_X%zeta_F,2) .or. &
            &n_r.ne.size(grid_X%zeta_F,3) .or. &
            &n_r.ne.size(grid_X%grp_r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and grp_r_F of grid_X need to have the &
                &correct dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(X_F,1) .or. n_zeta.ne.size(X_F,2) .or. &
            &n_r.ne.size(X_F,3) .or. n_t.ne.size(X_F,4)) then
            ierr = 1
            err_msg = 'X_F needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! set up multiplicative factor due to derivatives
        allocate(deriv_mult_factor(grp_n_r_X_loc))
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (abs(realpart(sqrt_X_val_norm)) .lt. imagpart(sqrt_X_val_norm)) then ! exploding, unstable
            if (imagpart(sqrt_X_val_norm).gt.0) &
                &sqrt_X_val_norm = - sqrt_X_val_norm                            ! exploding solution, not the decaying one
        end if
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)
        
        ! initialize X_F
        X_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,X%n_mod
                ! set up angular multiplicative factor for parallel derivative
                if (deriv_loc) then
                    if (use_pol_flux_X) then
                        deriv_mult_factor =  &
                            &iu*(X%n(jd)*eq%q_saf_FD(:,0)-X%m(jd))
                    else
                        deriv_mult_factor =  &
                            &iu*(X%n(jd)-X%m(jd)*eq%rot_t_FD(:,0))
                    end if
                else
                    deriv_mult_factor = 1._dp
                end if
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X_loc
                    ! add current mode
                    X_F(:,:,kd,id) = X_F(:,:,kd,id) + realpart(&
                        &exp(iu*(X%n(jd)*grid_X%zeta_F(:,:,kd)-&
                        &X%m(jd)*grid_X%theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor(kd) * X%vec(jd,kd,X_id))
                end do
            end do
        end do
    end function calc_real_X_arr
    integer function calc_real_X_ind(eq,grid_X,X,X_id,time,X_F,deriv) &
        &result(ierr)                                                           ! (time) individual version
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_real_X_ind'
        
        ! input / output
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: X_F(:,:,:)                                   ! normal component of perturbation
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        real(dp), allocatable :: X_F_arr(:,:,:,:)
        
        ! allocate array version of X_F
        allocate(X_F_arr(size(X_F,1),size(X_F,2),size(X_F,3),1))
        
        ! call array version
        ierr = calc_real_X_arr(eq,grid_X,X,X_id,[time],X_F_arr,deriv)
        CHCKERR('')
        
        ! copy array to individual X_F
        X_F = X_F_arr(:,:,:,1)
        
        ! deallocate array version of X_F
        deallocate(X_F_arr)
    end function calc_real_X_ind
    
    !! calculates the geodesic  component of the perturbation,  or optionally the
    !! parallel derivative.  The input is  given for a range  (r,theta,zeta)_F in
    !! Flux coordinates and for a range  in normalized time (1 corresponds to one
    !! period).
    !integer function calc_real_U_arr(X_vec,X_val,r_F,theta_F,zeta_F,time,U_F,&
        !&deriv) result(ierr)                                                    ! (time) array version
        !use X_vars, only: grp_n_r_X, size_X, n_X, m_X, U_X_0, U_X_1, DU_X_0, &
            !&DU_X_1, grp_r_X
        !use eq_vars, only: q_saf_FD, rot_t_FD
        !use utilities, only: calc_deriv
        !use num_vars, only: use_pol_flux_X
        
        !character(*), parameter :: rout_name = 'calc_real_U_arr'
        
        !! input / output
        !complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        !complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        !real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        !real(dp), intent(in) :: time(:)                                         ! time range
        !real(dp), intent(inout) :: U_F(:,:,:,:)                                 ! normal component of perturbation
        !logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        !! local variables
        !integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        !integer :: n_t                                                          ! number of time points
        !character(len=max_str_ln) :: err_msg                                    ! error message
        !integer :: id, jd, kd                                                   ! counter
        !complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        !logical :: deriv_loc                                                    ! local copy of deriv
        !complex(dp), allocatable :: deriv_mult_factor(:)                        ! multiplicative factor due to parallel derivative
        !complex(dp), allocatable :: DX_vec(:)                                   ! derivative of X_vec for a specific mode
        
        !! initialize ierr
        !ierr = 0
        
        !! set up array sizes and time points
        !n_theta = size(theta_F,1)
        !n_zeta = size(theta_F,2)
        !n_r = size(theta_F,3)
        !n_t = size(time)
        
        !! tests
        !if (n_theta.ne.size(zeta_F,1) .or. n_zeta.ne.size(zeta_F,2) .or. &
            !&n_r.ne.size(zeta_F,3) .or. n_r.ne.size(r_F)) then
            !ierr = 1
            !err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                !&dimensions'
            !CHCKERR(err_msg)
        !end if
        !if (n_theta.ne.size(U_F,1) .or. n_zeta.ne.size(U_F,2) .or. &
            !&n_r.ne.size(U_F,3) .or. n_t.ne.size(U_F,4)) then
            !ierr = 1
            !err_msg = 'U_F needs to have the correct dimensions'
            !CHCKERR(err_msg)
        !end if
        
        !! set up local copy of deriv
        !deriv_loc = .false.
        !if (present(deriv)) deriv_loc = deriv
        
        !! set up DX_vec
        !allocate(DX_vec(grp_n_r_X))
        
        !! set up multiplicative factor due to derivatives
        !allocate(deriv_mult_factor(grp_n_r_X))
        
        !! set up normalized sqrt(X_val)
        !sqrt_X_val_norm = sqrt(X_val)
        !if (abs(realpart(sqrt_X_val_norm)) .lt. imagpart(sqrt_X_val_norm)) then ! exploding, unstable
            !if (imagpart(sqrt_X_val_norm).gt.0) &
                !&sqrt_X_val_norm = - sqrt_X_val_norm                            ! exploding solution, not the decaying one
        !end if
        !sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)
        
        !! initialize U_F
        !U_F = 0._dp
        
        !! iterate over time steps
        !do id = 1,n_t
            !! iterate over all modes
            !do jd = 1,size_X
                !! set up angular multiplicative factor for parallel derivative
                !if (deriv_loc) then
                    !deriv_mult_factor = 1._dp
                !else
                    !deriv_mult_factor = 1._dp
                !end if
                
                !! set up normal derivative of X_vec
                !ierr = calc_deriv(X_vec(jd,:),DX_vec,grp_r_X(1:grp_n_r_X),1,2)
                !CHCKERR('')
                
                !! iterate over all normal points (of this group)
                !do kd = 1,grp_n_r_X
                    !! add current mode
                    !U_F(:,:,kd,id) = U_F(:,:,kd,id) + realpart(&
                        !&exp(iu*(n_X(jd)*zeta_F(:,:,kd)-&
                        !&m_X(jd)*theta_F(:,:,kd))) * &
                        !&exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        !&deriv_mult_factor(kd) * X_vec(jd,kd))
                !end do
                
                !! deallocate normal derivative of X_vec
                !deallocate(DX_vec)
            !end do
        !end do
    !end function calc_real_U_arr
    !integer function calc_real_U_ind(X_vec,X_val,r_F,theta_F,zeta_F,time,U_F,&
        !&deriv) result(ierr)                                                    ! (time) individual version
        
        !character(*), parameter :: rout_name = 'calc_real_U_ind'
        
        !! input / output
        !complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        !complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        !real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        !real(dp), intent(in) :: time                                            ! time range
        !real(dp), intent(inout) :: U_F(:,:,:)                                   ! normal component of perturbation
        !logical, intent(in), optional :: deriv                                  ! optional parallel derivative
        
        !! local variables
        !real(dp), allocatable :: U_F_arr(:,:,:,:)
        
        !! allocate array version of U_F
        !allocate(U_F_arr(size(U_F,1),size(U_F,2),size(U_F,3),1))
        
        !! call array version
        !ierr = calc_real_U_arr(X_vec,X_val,r_F,theta_F,zeta_F,[time],&
            !&U_F_arr,deriv)
        !CHCKERR('')
        
        !! copy array to individual U_F
        !U_F = U_F_arr(:,:,:,1)
        
        !! deallocate array version of U_F
        !deallocate(U_F_arr)
    !end function calc_real_U_ind
end module sol_ops
