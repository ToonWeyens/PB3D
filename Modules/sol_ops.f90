!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi

    implicit none
    private
    public calc_real_X, calc_real_U
    
    ! interfaces
    interface calc_real_X
        module procedure calc_real_X_arr, calc_real_X_ind
    end interface
    interface calc_real_U
        module procedure calc_real_U_arr, calc_real_U_ind
    end interface

contains
    ! calculates the normal  component of the perturbation,  or derivatives. The
    ! input is given for a range  (r,theta,zeta)_F in Flux coordinates and for a
    ! range in normalized time (1 corresponds to one period).
    ! Optionally, derivatives can be specified in [r,theta,z]
    integer function calc_real_X_arr(X_vec,X_val,r_F,theta_F,zeta_F,time,X_F,&
        &deriv) result(ierr)                                                    ! (time) array version
        use X_vars, only: grp_n_r_X, size_X, n_X, m_X
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: X_F(:,:,:,:)                                 ! normal component of perturbation
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        integer :: deriv_loc(3)                                                 ! local copy of deriv
        complex(dp) :: deriv_mult_factor                                        ! multiplicative factor due to angular derivatives
        real(dp), allocatable :: X_F_loc(:)                                     ! local version of X_F, used to calculate normal derivatives
        integer :: deriv_prec                                                   ! precision of normal derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(theta_F,1)
        n_zeta = size(theta_F,2)
        n_r = size(theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(zeta_F,1) .or. n_zeta.ne.size(zeta_F,2) .or. &
            &n_r.ne.size(zeta_F,3) .or. n_r.ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(X_F,1) .or. n_zeta.ne.size(X_F,2) .or. &
            &n_r.ne.size(X_F,3) .or. n_t.ne.size(X_F,4)) then
            ierr = 1
            err_msg = 'X_F needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(deriv)) then
            if (minval(deriv).lt.0) then
                ierr = 1
                err_msg = 'Only positive derivatives can be specified'
                CHCKERR(err_msg)
            end if
            if(deriv(1).gt.2) then
                ierr = 1
                err_msg = 'Normal derivatives cannot be of higher order than 2'
                CHCKERR(err_msg)
            end if
        end if
        
        ! set up local copy of deriv
        deriv_loc = [0,0,0]
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X_val)
        if (abs(realpart(sqrt_X_val_norm)) .lt. imagpart(sqrt_X_val_norm)) then ! exploding, unstable
            if (imagpart(sqrt_X_val_norm).gt.0) &
                &sqrt_X_val_norm = - sqrt_X_val_norm                            ! exploding solution, not the decaying one
        end if
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)
        
        ! if normal derivatives, allocate helper variable
        if (deriv_loc(1).gt.0) then
            allocate(X_F_loc(n_r))
            if (deriv_loc(1).eq.1) deriv_prec = 2
            if (deriv_loc(1).eq.2) deriv_prec = 1
        end if
        
        ! initialize X_F
        X_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,size_X
                ! set up angular multiplicative factor for derivatives
                deriv_mult_factor =  &
                    &(-iu*m_X(jd))**deriv_loc(2) * (iu*n_X(jd))**deriv_loc(3)
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X
                    ! add current mode
                    X_F(:,:,kd,id) = X_F(:,:,kd,id) + realpart(&
                        &exp(iu*(n_X(jd)*zeta_F(:,:,kd)-&
                        &m_X(jd)*theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor * X_vec(jd,kd))
                end do
            end do
            
            ! perform normal derivatives on global solution
            if (deriv_loc(1).gt.0) then
                do jd = 1,n_zeta
                    do kd = 1,n_theta
                        X_F_loc = X_F(kd,jd,:,id)
                        ierr = calc_deriv(X_F_loc,X_F(kd,jd,:,id),r_F,&
                            &deriv_loc(1),deriv_prec)
                        CHCKERR('')
                    end do
                end do
            end if
        end do
        
        ! deallocate helper variables
        if (deriv_loc(1).gt.0) then
            deallocate(X_F_loc)
        end if
    end function calc_real_X_arr
    integer function calc_real_X_ind(X_vec,X_val,r_F,theta_F,zeta_F,time,X_F,&
        &deriv) result(ierr)                                                    ! (time) individual version
        
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(in) :: time                                            ! time range
        real(dp), intent(inout) :: X_F(:,:,:)                                   ! normal component of perturbation
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        real(dp), allocatable :: X_F_arr(:,:,:,:)
        
        ! allocate array version of X_F
        allocate(X_F_arr(size(X_F,1),size(X_F,2),size(X_F,3),1))
        
        ! call array version
        ierr = calc_real_X_arr(X_vec,X_val,r_F,theta_F,zeta_F,[time],&
            &X_F_arr,deriv)
        CHCKERR('')
        
        ! copy array to individual X_F
        X_F = X_F_arr(:,:,:,1)
        
        ! deallocate array version of X_F
        deallocate(X_F_arr)
    end function calc_real_X_ind
    
    ! calculates the geodesic component of the perturbation, or derivatives. The
    ! input is given for a range  (r,theta,zeta)_F in Flux coordinates and for a
    ! range in normalized time (1 corresponds to one period).
    ! Optionally, derivatives can be specified in [r,theta,z]
    integer function calc_real_U_arr(X_vec,X_val,r_F,theta_F,zeta_F,time,U_F,&
        &deriv) result(ierr)                                                    ! (time) array version
        use X_vars, only: grp_n_r_X, size_X, n_X, m_X
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: U_F(:,:,:,:)                                 ! normal component of perturbation
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        integer :: deriv_loc(3)                                                 ! local copy of deriv
        complex(dp) :: deriv_mult_factor                                        ! multiplicative factor due to angular derivatives
        real(dp), allocatable :: X_F_loc(:)                                     ! local version of X_F, used to calculate normal derivatives
        integer :: deriv_prec                                                   ! precision of normal derivatives
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(theta_F,1)
        n_zeta = size(theta_F,2)
        n_r = size(theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(zeta_F,1) .or. n_zeta.ne.size(zeta_F,2) .or. &
            &n_r.ne.size(zeta_F,3) .or. n_r.ne.size(r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and r_F need to have the correct &
                &dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(U_F,1) .or. n_zeta.ne.size(U_F,2) .or. &
            &n_r.ne.size(U_F,3) .or. n_t.ne.size(U_F,4)) then
            ierr = 1
            err_msg = 'U_F needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if (present(deriv)) then
            if (minval(deriv).lt.0) then
                ierr = 1
                err_msg = 'Only positive derivatives can be specified'
                CHCKERR(err_msg)
            end if
            if(deriv(1).gt.2) then
                ierr = 1
                err_msg = 'Normal derivatives cannot be of higher order than 2'
                CHCKERR(err_msg)
            end if
        end if
        
        ! set up local copy of deriv
        deriv_loc = [0,0,0]
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X_val)
        if (abs(realpart(sqrt_X_val_norm)) .lt. imagpart(sqrt_X_val_norm)) then ! exploding, unstable
            if (imagpart(sqrt_X_val_norm).gt.0) &
                &sqrt_X_val_norm = - sqrt_X_val_norm                            ! exploding solution, not the decaying one
        end if
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)
        
        ! if normal derivatives, allocate helper variable
        if (deriv_loc(1).gt.0) then
            allocate(X_F_loc(n_r))
            if (deriv_loc(1).eq.1) deriv_prec = 2
            if (deriv_loc(1).eq.2) deriv_prec = 1
        end if
        
        ! initialize X_F
        U_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,size_X
                ! set up angular multiplicative factor for derivatives
                deriv_mult_factor =  &
                    &(-iu*m_X(jd))**deriv_loc(2) * (iu*n_X(jd))**deriv_loc(3)
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X
                    ! add current mode
                    U_F(:,:,kd,id) = U_F(:,:,kd,id) + realpart(&
                        &exp(iu*(n_X(jd)*zeta_F(:,:,kd)-&
                        &m_X(jd)*theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor * X_vec(jd,kd))
                end do
            end do
            
            ! perform normal derivatives on global solution
            if (deriv_loc(1).gt.0) then
                do jd = 1,n_zeta
                    do kd = 1,n_theta
                        X_F_loc = U_F(kd,jd,:,id)
                        ierr = calc_deriv(X_F_loc,U_F(kd,jd,:,id),r_F,&
                            &deriv_loc(1),deriv_prec)
                        CHCKERR('')
                    end do
                end do
            end if
        end do
        
        ! deallocate helper variables
        if (deriv_loc(1).gt.0) then
            deallocate(X_F_loc)
        end if
    end function calc_real_U_arr
    integer function calc_real_U_ind(X_vec,X_val,r_F,theta_F,zeta_F,time,U_F,&
        &deriv) result(ierr)                                                    ! (time) individual version
        
        character(*), parameter :: rout_name = 'calc_real_X'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! MPI Eigenvector
        complex(dp), intent(in) :: X_val                                        ! Eigenvalue
        real(dp), intent(in) :: r_F(:), theta_F(:,:,:), zeta_F(:,:,:)           ! Flux (perturbation) coords.
        real(dp), intent(in) :: time                                            ! time range
        real(dp), intent(inout) :: U_F(:,:,:)                                   ! normal component of perturbation
        integer, intent(in), optional :: deriv(3)                               ! optional derivatives in angular coordinates
        
        ! local variables
        real(dp), allocatable :: X_F_arr(:,:,:,:)
        
        ! allocate array version of U_F
        allocate(X_F_arr(size(U_F,1),size(U_F,2),size(U_F,3),1))
        
        ! call array version
        ierr = calc_real_X_arr(X_vec,X_val,r_F,theta_F,zeta_F,[time],&
            &X_F_arr,deriv)
        CHCKERR('')
        
        ! copy array to individual X_F
        U_F = X_F_arr(:,:,:,1)
        
        ! deallocate array version of X_F
        deallocate(X_F_arr)
    end function calc_real_U_ind
end module sol_ops
