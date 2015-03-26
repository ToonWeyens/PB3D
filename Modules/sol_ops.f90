!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, max_str_ln, pi
    use output_ops, only: print_GP_2D

    implicit none
    private
    public calc_real_XUQ, calc_real_U, calc_real_Qn, calc_real_Qg
    
    ! interfaces
    interface calc_real_XUQ
        module procedure calc_real_XUQ_arr, calc_real_XUQ_ind
    end interface
    interface calc_real_U
        module procedure calc_real_U_arr, calc_real_U_ind
    end interface
    interface calc_real_Qn
        module procedure calc_real_Qn_arr, calc_real_Qn_ind
    end interface
    interface calc_real_Qg
        module procedure calc_real_Qg_arr, calc_real_Qg_ind
    end interface

contains
    ! calculates the  normal or geodesic  component, or optionally  the parallel
    ! derivative of the plasma perturbation  or the normal or geodesic component
    ! of the magnetic field perturbation in the perturbation grid for a range in
    ! normalized  time  (1  corresponds  to  one  period).  The  variable  style
    ! determines which one:
    !   - style = 1: X (supports parallel derivative)
    !   - style = 2: U (supports parallel derivative)
    !   - style = 3: Qn
    !   - style = 4: Qg
    ! For Qn  and Qg,  the metric  variables have  to be  provided as  well. The
    ! output is  given in the  output grid for the  normal part, onto  which the
    ! normal part  of the equilibrium  and perturbation grids  are interpolated.
    ! The angular part of the output is given by the equilibrium grid.
    ! Note: The angular part of the output grid is neglected as it is assumed to
    ! be a 1D grid.
    integer function calc_real_XUQ_arr(grid_eq,eq,grid_X,X,grid_XUQ,X_id,style,&
        &time,XUQ,met,deriv) result(ierr)                                       ! (time) array version
        use grid_vars, only: grid_type
        use eq_vars, only: eq_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        use num_vars, only: use_pol_flux_F, grp_rank, grp_n_procs
        use utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'calc_real_XUQ_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        type(grid_type), intent(in) :: grid_XUQ                                 ! output grid
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: style                                            ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: XUQ(:,:,:,:)                                 ! normal component of perturbation
        type(metric_type), optional, intent(in) :: met                          ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        complex(dp), allocatable :: deriv_mult_factor(:)                        ! multiplicative factor due to parallel derivative
        complex(dp), pointer :: XU_0(:,:,:), XU_1(:,:,:)                        ! either 1 and 0, U_0 and U_1 or DU_0 and DU_1
        real(dp), allocatable :: grp_r_X(:)                                     ! grp_r_F of XUQ interpolated in X grid
        real(dp), allocatable :: grp_r_eq(:)                                    ! grp_r_F of XUQ interpolated in eq grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_t
        n_t = size(time)
        
        ! tests
        if (size(grid_XUQ%theta_F,1).ne.size(XUQ,1) .or. &
            &size(grid_XUQ%theta_F,2).ne.size(XUQ,2) .or. &
            &size(grid_XUQ%theta_F,3).ne.size(XUQ,3) .or. &
            &n_t.ne.size(XUQ,4)) then
            ierr = 1
            err_msg = 'XUQ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if ((style.eq.3 .or. style.eq.4).and..not.present(met)) then
            ierr = 1
            err_msg = 'When calculating Q, also metric variables have to be &
                &provided'
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
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! calculate  normal  interpolation  tables for  equilibrium  grid  (this
        ! concerns eq variables and met variables)
        allocate(grp_r_eq(grid_XUQ%grp_n_r))
        do kd = 1,grid_XUQ%grp_n_r
            call con2dis(grid_XUQ%grp_r_F(kd),grp_r_eq(kd),grid_eq%grp_r_F)
        end do
        
        ! calculate  normal interpolation  tables  for  perturbation grid  (this
        ! concerns eq variables and met variables)
        allocate(grp_r_X(grid_XUQ%grp_n_r))
        do kd = 1,grid_XUQ%grp_n_r
            call con2dis(grid_XUQ%grp_r_F(kd),grp_r_X(kd),grid_X%grp_r_F)
        end do
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,X%n_mod
                ! set up angular multiplicative factor for parallel derivative
                write(*,*) '!!!!!!!!! NEED TO INTERPOLATE EQ GRID !!!!!!!!!!'
                if (deriv_loc) then
                    if (use_pol_flux_F) then
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
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + realpart(&
                        &exp(iu*(X%n(jd)*grid_X%zeta_F(:,:,kd)-&
                        &X%m(jd)*grid_X%theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor(kd) * X%vec(jd,kd,X_id))
                end do
            end do
        end do
    end function calc_real_XUQ_arr
    integer function calc_real_XUQ_ind(grid_eq,eq,grid_X,X,grid_XUQ,X_id,&
        &style,time,XUQ,met,deriv) result(ierr)                                                ! (time) individual version
        use grid_vars, only: grid_type
        use eq_vars, only: eq_type
        use metric_vars, only: metric_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_real_XUQ_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibirum grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        type(grid_type), intent(in) :: grid_XUQ                                 ! output grid
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: style                                            ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: XUQ(:,:,:)                                   ! normal component of perturbation
        type(metric_type), optional, intent(in) :: met                          ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        real(dp), allocatable :: XUQ_arr(:,:,:,:)
        
        ! allocate array version of XUQ
        allocate(XUQ_arr(size(XUQ,1),size(XUQ,2),size(XUQ,3),1))
        
        ! call array version
        ierr = calc_real_XUQ_arr(grid_eq,eq,grid_X,X,grid_XUQ,X_id,style,&
            &[time],XUQ_arr,met,deriv)
        CHCKERR('')
        
        ! copy array to individual XUQ
        XUQ = XUQ_arr(:,:,:,1)
        
        ! deallocate array version of XUQ
        deallocate(XUQ_arr)
    end function calc_real_XUQ_ind
    
    ! calculates the geodesic  component of the perturbation,  or optionally the
    ! parallel derivative for  a range in normalized time (1  corresponds to one
    ! period).
    integer function calc_real_U_arr(grid_eq,eq,grid_X,X,X_id,time,U_F,deriv) &
        &result(ierr)                                                           ! (time) array version
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        use num_vars, only: use_pol_flux_F, grp_rank, grp_n_procs
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'calc_real_U_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: U_F(:,:,:,:)                                 ! geodesic component of perturbation
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        complex(dp), allocatable :: deriv_mult_factor(:)                        ! multiplicative factor due to parallel derivative
        complex(dp), allocatable :: DX_vec(:)                                   ! derivative of X_vec for a specific mode
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        complex(dp), pointer :: U_0(:,:,:), U_1(:,:,:)                          ! either U_0 and U_1 or DU_0 and DU_1
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(grid_eq%theta_F,1)
        n_zeta = size(grid_eq%theta_F,2)
        n_r = size(grid_eq%theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(grid_eq%zeta_F,1) .or. &
            &n_zeta.ne.size(grid_eq%zeta_F,2) .or. &
            &n_r.ne.size(grid_eq%zeta_F,3) .or. &
            &n_r.ne.size(grid_eq%grp_r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and grp_r_F of grid_eq need to have the &
                &correct dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(U_F,1) .or. n_zeta.ne.size(U_F,2) .or. &
            &n_r.ne.size(U_F,3) .or. n_t.ne.size(U_F,4)) then
            ierr = 1
            err_msg = 'U_F needs to have the correct dimensions'
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
        
        ! set up DX_vec
        allocate(DX_vec(grp_n_r_X_loc))
        
        ! set up multiplicative factor due to derivatives
        allocate(deriv_mult_factor(grp_n_r_X_loc))
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! initialize U_F
        U_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,X%n_mod
                ! set up angular multiplicative factor for parallel derivative
                if (deriv_loc) then
                    if (use_pol_flux_F) then
                        deriv_mult_factor =  &
                            &iu*(X%n(jd)*eq%q_saf_FD(:,0)-X%m(jd))
                    else
                        deriv_mult_factor =  &
                            &iu*(X%n(jd)-X%m(jd)*eq%rot_t_FD(:,0))
                    end if
                    U_0 => X%U_0(:,:,:,jd)
                    U_1 => X%U_1(:,:,:,jd)
                else
                    deriv_mult_factor = 1._dp
                    U_0 => X%DU_0(:,:,:,jd)
                    U_1 => X%DU_1(:,:,:,jd)
                end if
                
                ! set up normal derivative of X_vec
                ierr = calc_deriv(X%vec(jd,:,X_id),DX_vec,&
                    &grid_X%grp_r_F(1:grp_n_r_X_loc),1,2)
                CHCKERR('')
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X_loc
                    ! add current mode
                    U_F(:,:,kd,id) = U_F(:,:,kd,id) + realpart(&
                        &exp(iu*(X%n(jd)*grid_X%zeta_F(:,:,kd)-&
                        &X%m(jd)*grid_X%theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor(kd) * &
                        &(U_0(:,:,kd)*X%vec(jd,kd,X_id) + &
                        &U_1(:,:,kd)*DX_vec(kd)))
                end do
            end do
        end do
    end function calc_real_U_arr
    integer function calc_real_U_ind(grid_eq,eq,grid_X,X,X_id,time,U_F,deriv) &
        &result(ierr)                                                           ! (time) individual version
        use eq_vars, only: eq_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_real_U_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! perturbation grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: U_F(:,:,:)                                   ! geodesic component of perturbation
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        real(dp), allocatable :: U_F_arr(:,:,:,:)
        
        ! allocate array version of U_F
        allocate(U_F_arr(size(U_F,1),size(U_F,2),size(U_F,3),1))
        
        ! call array version
        ierr = calc_real_U_arr(grid_eq,eq,grid_X,X,X_id,[time],U_F_arr,deriv)
        CHCKERR('')
        
        ! copy array to individual U_F
        U_F = U_F_arr(:,:,:,1)
        
        ! deallocate array version of U_F
        deallocate(U_F_arr)
    end function calc_real_U_ind
    
    ! calculates   the  normal   component  of   the  modified   magnetic  field
    ! perturbation  Q. The  input is  given for  a range  in normalized  time (1
    ! corresponds to one period).
    integer function calc_real_Qn_arr(grid_eq,eq,met,grid_X,X,X_id,time,Qn_F) &
        &result(ierr)  ! (time) array version
        use eq_vars, only: eq_type
        use metric_vars, only: metric_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        use num_vars, only: use_pol_flux_F, grp_rank, grp_n_procs
        
        character(*), parameter :: rout_name = 'calc_real_Qn_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: Qn_F(:,:,:,:)                                ! normal component of Q
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        complex(dp), allocatable :: deriv_mult_factor(:)                        ! multiplicative factor due to parallel derivative
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(grid_eq%theta_F,1)
        n_zeta = size(grid_eq%theta_F,2)
        n_r = size(grid_eq%theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(grid_eq%zeta_F,1) .or. &
            &n_zeta.ne.size(grid_eq%zeta_F,2) .or. &
            &n_r.ne.size(grid_eq%zeta_F,3) .or. &
            &n_r.ne.size(grid_eq%grp_r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and grp_r_F of grid_eq need to have the &
                &correct dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(Qn_F,1) .or. n_zeta.ne.size(Qn_F,2) .or. &
            &n_r.ne.size(Qn_F,3) .or. n_t.ne.size(Qn_F,4)) then
            ierr = 1
            err_msg = 'Qn_F needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! set up multiplicative factor due to derivatives
        allocate(deriv_mult_factor(grp_n_r_X_loc))
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! initialize Qn_F
        Qn_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,X%n_mod
                ! set up angular multiplicative factor for parallel derivative
                if (use_pol_flux_F) then
                    deriv_mult_factor =  &
                        &iu*(X%n(jd)*eq%q_saf_FD(:,0)-X%m(jd))
                else
                    deriv_mult_factor =  &
                        &iu*(X%n(jd)-X%m(jd)*eq%rot_t_FD(:,0))
                end if
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X_loc
                    ! add current mode
                    Qn_F(:,:,kd,id) = Qn_F(:,:,kd,id) + realpart(&
                        &exp(iu*(X%n(jd)*grid_X%zeta_F(:,:,kd)-&
                        &X%m(jd)*grid_X%theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &deriv_mult_factor(kd) * X%vec(jd,kd,X_id))
                end do
            end do
            
            ! divide by the jacobian
            Qn_F(:,:,:,id) = Qn_F(:,:,:,id)/met%jac_FD(:,:,:,0,0,0)
        end do
    end function calc_real_Qn_arr
    integer function calc_real_Qn_ind(grid_eq,eq,met,grid_X,X,X_id,time,Qn_F) &
        &result(ierr)                                                           ! (time) individual version
        use eq_vars, only: eq_type
        use metric_vars, only: metric_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_real_Qn_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: Qn_F(:,:,:)                                  ! normal component of Q
        
        ! local variables
        real(dp), allocatable :: Qn_F_arr(:,:,:,:)
        
        ! allocate array version of Qn_F
        allocate(Qn_F_arr(size(Qn_F,1),size(Qn_F,2),size(Qn_F,3),1))
        
        ! call array version
        ierr = calc_real_Qn_arr(grid_eq,eq,met,grid_X,X,X_id,[time],Qn_F_arr)
        CHCKERR('')
        
        ! copy array to individual Qn_F
        Qn_F = Qn_F_arr(:,:,:,1)
        
        ! deallocate array version of Qn_F
        deallocate(Qn_F_arr)
    end function calc_real_Qn_ind
    
    ! calculates  the   geodesic  component  of  the   modified  magnetic  field
    ! perturbation Q. The input is given for a range in normalized time (1
    ! corresponds to one period).
    integer function calc_real_Qg_arr(grid_eq,met,grid_X,X,X_id,time,Qg_F) &
        &result(ierr)                                                           ! (time) array version
        use metric_vars, only: metric_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        use num_vars, only: grp_rank, grp_n_procs
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'calc_real_Qg_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: Qg_F(:,:,:,:)                                ! geodesic component of Q
        
        ! local variables
        integer :: n_r, n_theta, n_zeta                                         ! dimensions of the grid
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        complex(dp), allocatable :: DX_vec(:)                                   ! derivative of X_vec for a specific mode
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up array sizes and time points
        n_theta = size(grid_eq%theta_F,1)
        n_zeta = size(grid_eq%theta_F,2)
        n_r = size(grid_eq%theta_F,3)
        n_t = size(time)
        
        ! tests
        if (n_theta.ne.size(grid_eq%zeta_F,1) .or. &
            &n_zeta.ne.size(grid_eq%zeta_F,2) .or. &
            &n_r.ne.size(grid_eq%zeta_F,3) .or. &
            &n_r.ne.size(grid_eq%grp_r_F)) then
            ierr = 1
            err_msg = 'theta_F, zeta_F and grp_r_F of grid_eq need to have the &
                &correct dimensions'
            CHCKERR(err_msg)
        end if
        if (n_theta.ne.size(Qg_F,1) .or. n_zeta.ne.size(Qg_F,2) .or. &
            &n_r.ne.size(Qg_F,3) .or. n_t.ne.size(Qg_F,4)) then
            ierr = 1
            err_msg = 'Qg_F needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! set up DX_vec
        allocate(DX_vec(grp_n_r_X_loc))
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! initialize Qg_F
        Qg_F = 0._dp
        
        ! iterate over time steps
        do id = 1,n_t
            ! iterate over all modes
            do jd = 1,X%n_mod
                ! set up normal derivative of X_vec
                ierr = calc_deriv(X%vec(jd,:,X_id),DX_vec,&
                    &grid_X%grp_r_F(1:grp_n_r_X_loc),1,2)
                CHCKERR('')
                
                ! iterate over all normal points (of this group)
                do kd = 1,grp_n_r_X_loc
                    ! add current mode
                    Qg_F(:,:,kd,id) = Qg_F(:,:,kd,id) + realpart(&
                        &exp(iu*(X%n(jd)*grid_X%zeta_F(:,:,kd)-&
                        &X%m(jd)*grid_X%theta_F(:,:,kd))) * &
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * &
                        &((X%DU_0(:,:,kd,jd)-X%extra1(:,:,kd))*&
                        &X%vec(jd,kd,X_id) + X%DU_1(:,:,kd,jd)*DX_vec(kd)))
                end do
            end do
            
            ! divide by the jacobian
            Qg_F(:,:,:,id) = Qg_F(:,:,:,id)/met%jac_FD(:,:,:,0,0,0)
        end do
    end function calc_real_Qg_arr
    integer function calc_real_Qg_ind(grid_eq,met,grid_X,X,X_id,time,Qg_F) &
        &result(ierr)                                                           ! (time) individual version
        use metric_vars, only: metric_type
        use grid_vars, only: grid_type
        use X_vars, only: X_type
        
        character(*), parameter :: rout_name = 'calc_real_Qg_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(metric_type), intent(in) :: met                                    ! metric variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: Qg_F(:,:,:)                                  ! geodesic component of Q
        
        ! local variables
        real(dp), allocatable :: Qg_F_arr(:,:,:,:)
        
        ! allocate array version of Qg_F
        allocate(Qg_F_arr(size(Qg_F,1),size(Qg_F,2),size(Qg_F,3),1))
        
        ! call array version
        ierr = calc_real_Qg_arr(grid_eq,met,grid_X,X,X_id,[time],Qg_F_arr)
        CHCKERR('')
        
        ! copy array to individual Qg_F
        Qg_F = Qg_F_arr(:,:,:,1)
        
        ! deallocate array version of Qg_F
        deallocate(Qg_F_arr)
    end function calc_real_Qg_ind
end module sol_ops
