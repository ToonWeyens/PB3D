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
    public prepare_X, solve_EV_system, plot_X_vec

contains
    ! prepare the matrix elements by calculating KV and PV, which then will have
    ! to be integrated, with a complex exponential weighting function
    subroutine prepare_X
        use X_vars, only: init_X, calc_PV, calc_KV, calc_U, calc_extra, &
            &calc_V_int, dealloc_X_vars, &
            &PV0, PV1, PV2, KV0, KV1, KV2, PV_int, KV_int
        
        call writo('Calculating table of field line averages &
            &<V e^i(k-m)theta>')
        
        call lvl_ud(1)
        
        ! initialize the variables
        call writo('Initalizing variables...')
        call lvl_ud(1)
        call init_X
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
        
        ! Calculate PV_int = <PVi e^(k-m)theta>
        call writo('Taking field average of PV')
        call lvl_ud(1)
        call calc_V_int(PV0,PV_int(:,:,:,1))
        call calc_V_int(PV1,PV_int(:,:,:,2))
        call calc_V_int(PV2,PV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! Calculate KV_int = <KVi e^(k-m)theta>
        call writo('Taking field average of KV')
        call lvl_ud(1)
        call calc_V_int(KV0,KV_int(:,:,:,1))
        call calc_V_int(KV1,KV_int(:,:,:,2))
        call calc_V_int(KV2,KV_int(:,:,:,3))
        call lvl_ud(-1)
        
        ! deallocate equilibrium variables
        call writo('deallocating unused variables')
        call dealloc_X_vars
        
        call lvl_ud(-1)
        
        call writo('Done calculating')
    end subroutine prepare_X
    
    ! set-up and  solve the  EV system  by discretizing  the equations  in n_r_X
    ! normal points,  making use of  PV0, PV1 and  PV2, interpolated in  the n_r
    ! (equilibrium) values
    integer function solve_EV_system(A_terms,B_terms,A_info,B_info) result(ierr)
        use num_vars, only: EV_style
        use str_ops, only: i2str
        use slepc_ops, only: solve_EV_system_slepc
        
        character(*), parameter :: rout_name = 'solve_EV_system'
        
        ! input / output
        complex(dp), intent(inout), allocatable, optional :: A_terms(:,:,:,:)   ! termj_int of matrix A from previous Richardson loop
        complex(dp), intent(inout), allocatable, optional :: B_terms(:,:,:,:)   ! termj_int of matrix B from previous Richardson loop
        real(dp), intent(inout), optional :: A_info(2)                          ! info about A_terms: min of r and step_size
        real(dp), intent(inout), optional :: B_info(2)                          ! info about B_terms: min of r and step_size
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        select case (EV_style)
            case(1)                                                             ! slepc solver for EV problem
                ! solve the system
                ierr = solve_EV_system_slepc(A_terms,B_terms,A_info,B_info)
                CHCKERR('')
                
            case default
                err_msg = 'No EV solver style associated with '//&
                    &trim(i2str(EV_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function solve_EV_system
    
    ! plots an eigenfunction
    ! [MPI] Collective call
    integer function plot_X_vec(X_vec) result(ierr)
        use num_vars, only: grp_rank, min_r_X, max_r_X
        use X_vars, only: n_r_X, n_m_X
        use MPI_ops, only: get_ser_X_vec
        
        character(*), parameter :: rout_name = 'plot_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(n_m_X,n_r_X)                           ! MPI Eigenvector
        
        ! local variables
        real(dp), allocatable :: x_plot(:,:)                                    ! x_axis of plot
        complex(dp), allocatable :: ser_X_vec(:,:)                              ! serial version of X_vec
        integer :: id                                                           ! counter
        integer :: jd, kd                                                       ! counters
        real(dp), allocatable :: max_of_modes(:)                                ! maximum of each mode
        real(dp) :: current_magn                                                ! maximum of each mode
        real(dp), allocatable :: max_of_modes_r(:)                              ! flux surface where max of mode occurs
        
        ! initialize ierr
        ierr = 0
        
        ! initialize some things
        if (grp_rank.eq.0) then                                                 ! only for group masters
            allocate(ser_X_vec(n_m_X,n_r_X))
            allocate(max_of_modes(n_m_X))
            allocate(max_of_modes_r(n_m_X))
            max_of_modes = 0.0_dp
            max_of_modes_r = 0.0_dp
        else
            allocate(ser_X_vec(0,0))
        end if
        
        ! get serial version of X_vec on group master
        ierr = get_ser_X_vec(X_vec,ser_X_vec,n_m_X)
        CHCKERR('')
        
        ! calculations for group master
        if (grp_rank.eq.0) then                                                 ! only for group masters
            ! loop over all normal points of all modes in perturbation grid
            do kd = 1,n_r_X
                do jd = 1,n_m_X
                    ! check for maximum of mode jd and normal point jd
                    current_magn = sqrt(realpart(ser_X_vec(jd,kd))**2&
                        &+imagpart(ser_X_vec(jd,kd))**2)
                    if (current_magn.gt.max_of_modes(jd)) then
                        max_of_modes(jd) = current_magn
                        max_of_modes_r(jd) = min_r_X + (max_r_X-min_r_X) &
                            &* (kd-1.0)/(n_r_X-1.0)
                    end if
                end do
            end do
            call print_GP_2D('maximum of the modes','',max_of_modes)
            call print_GP_2D('place of maximum of the modes','',&
                &max_of_modes_r)
            
            deallocate(max_of_modes)
            deallocate(max_of_modes_r)
            
            ! set up x-axis
            allocate(x_plot(1:n_m_X,1:n_r_X))
            do id = 1,n_r_X
                x_plot(:,id) = min_r_X + (id-1.0)/(n_r_X-1)*(max_r_X-min_r_X)
            end do
            call print_GP_2D('norm of solution','',&
                &transpose(sqrt(realpart(ser_X_vec(:,:))**2 + &
                &imagpart(ser_X_vec(:,:))**2)),x=transpose(x_plot))
            
            ! deallocate variables
            deallocate(ser_X_vec)
            deallocate(x_plot)
        end if
    end function plot_X_vec
end module
