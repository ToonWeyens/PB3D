!------------------------------------------------------------------------------!
!   Calculate the matrix  elements due to the Plasma of  the system of equations
!   that has to be solved as described in [ADD REF]
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
    use num_vars, only: dp, iu, mu_0, max_str_ln
    use output_ops, only: lvl_ud, writo, print_ar_2, print_GP_2D
    use str_ops, only: i2str, r2strt

    implicit none
    private
    public prepare_matrix_X, solve_EV_system

contains
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    integer function prepare_matrix_X() result(ierr)
        use X_vars, only: init_X_ops, calc_PV, calc_KV, calc_U, calc_extra
        
        character(*), parameter :: rout_name = 'prepare_matrix_X'
        
        ! initialize ierr
        ierr = 0
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        ierr = init_X_ops()
        CHCKERR('')
        call lvl_ud(-1)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        call calc_U
        call lvl_ud(-1)
        
        ! calculate extra equilibrium quantities
        call writo('Calculating extra equilibrium quantities...')
        call lvl_ud(1)
        call calc_extra
        call lvl_ud(-1)
        
        ! Calculate PV0, PV1  and PV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV0, PV1 and PV2...')
        call lvl_ud(1)
        call calc_PV
        call lvl_ud(-1)
        
        ! Calculate KV0, KV1  and KV2 for all (k,m) pairs  and n_r (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV0, KV1 and KV2...')
        call lvl_ud(1)
        call calc_KV
        call lvl_ud(-1)
    end function prepare_matrix_X
    
    ! set-up and  solve the  EV system  by discretizing  the equations  in n_r_X
    ! normal points,  making use of  PV0, PV1 and  PV2, interpolated in  the n_r
    ! (equilibrium) values
    integer function solve_EV_system() result(ierr)
        use X_vars, only: n_r_X, X_vec, X_val, min_m_X, max_m_X, min_r, max_r
        use num_vars, only: EV_style, group_rank
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: x_plot(:,:)                                    ! x_axis of plot
        integer :: id, jd, kd                                                   ! counters
        real(dp), allocatable :: max_of_modes(:)                                ! maximum of each mode
        real(dp) :: current_magn                                                ! maximum of each mode
        real(dp), allocatable :: max_of_modes_r(:)                              ! flux surface where max of mode occurs
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                ierr = solve_EV_system_slepc()
                CHCKERR('')
                
                if (group_rank.eq.0) then
                    call writo('plotting results for first mode with eigenvalue '//&
                        &trim(r2strt(realpart(X_val(1))))//' + '//&
                        &(r2strt(imagpart(X_val(1)))))
                    call lvl_ud(1)
                
                    ! initialize max. of modes
                    allocate(max_of_modes(min_m_X:max_m_X))
                    allocate(max_of_modes_r(min_m_X:max_m_X))
                    max_of_modes = 0.0_dp
                    max_of_modes_r = 0.0_dp
                            
                    ! loop over all normal points of all modes in perturbation grid
                    do kd = 1,n_r_X
                        do jd = min_m_X,max_m_X
                            ! check for maximum of mode jd and normal point jd
                            current_magn = sqrt(realpart(X_vec(jd,kd,1))**2&
                                &+imagpart(X_vec(jd,kd,1))**2)
                            if (current_magn.gt.max_of_modes(jd)) then
                                max_of_modes(jd) = current_magn
                                max_of_modes_r(jd) = min_r + (max_r-min_r) &
                                    &* (kd-1.0)/(n_r_X-1.0)
                            end if
                        end do
                    end do
                    !call print_GP_2D('maximum of the modes','',max_of_modes)
                    !call print_GP_2D('place of maximum of the modes','',&
                        !&max_of_modes_r)
                    
                    deallocate(max_of_modes)
                    deallocate(max_of_modes_r)
                    
                    ! set up x-axis
                    allocate(x_plot(1:n_r_X,min_m_X:max_m_X))
                    do id = 1,n_r_X
                        x_plot(id,:) = min_r + (id-1.0)/(n_r_X-1)*(max_r-min_r)
                    end do
                    !call print_GP_2D('norm of solution','',&
                        !&transpose(sqrt(realpart(X_vec(:,:,1))**2 + &
                        !&imagpart(X_vec(:,:,1))**2)),x=x_plot)
                    deallocate(x_plot)
                    
                    call lvl_ud(-1)
                    
                end if
                
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
end module
