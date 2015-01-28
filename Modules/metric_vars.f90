!------------------------------------------------------------------------------!
!   Variables that have to do with the metric elements                         !
!------------------------------------------------------------------------------!
module metric_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, max_deriv, max_str_ln
    use str_ops, only: r2str, i2str
    
    implicit none
    private
    public dealloc_metric, dealloc_metric_final, &
        &jac_F, h_F, g_F, g_C, g_FD, h_FD, jac_E, jac_FD, T_FE, T_EF, &
        &det_T_FE, det_T_EF, h_E, g_E, T_VC, det_T_VC, jac_C

    ! upper (h) and lower (g) metric factors
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                                 ! in the C(ylindrical) coordinate system
    real(dp), allocatable :: g_E(:,:,:,:,:,:,:), h_E(:,:,:,:,:,:,:)             ! in the E(quilibrium) coordinate system
    real(dp), allocatable :: g_F(:,:,:,:,:,:,:), h_F(:,:,:,:,:,:,:)             ! in the F(lux) coordinate system with derivatves in the V(MEC) system
    real(dp), allocatable :: g_FD(:,:,:,:,:,:,:), h_FD(:,:,:,:,:,:,:)           ! in the F(lux) coordinate system with derivatives in the F(lux) system
    ! upper and lower transformation matrices
    ! (index 1,2: along B, perp to flux surfaces, index 4,5: 3x3 matrix)
    real(dp), allocatable :: T_VC(:,:,:,:,:,:,:)                                ! C(ylindrical) to V(MEC) (lower)
    real(dp), allocatable :: T_EF(:,:,:,:,:,:,:)                                ! E(quilibrium) to F(lux) (upper)
    real(dp), allocatable, target :: T_FE(:,:,:,:,:,:,:)                        ! E(quilibrium) to F(lux) (lower)
    real(dp), allocatable :: det_T_VC(:,:,:,:,:)                                ! determinant of T_VC
    real(dp), allocatable :: det_T_EF(:,:,:,:,:)                                ! determinant of T_EF
    real(dp), allocatable :: det_T_FE(:,:,:,:,:)                                ! determinant of T_FE
    real(dp), allocatable :: jac_C(:,:,:,:,:)                                   ! jacobian of C(ylindrical) coordinate system
    real(dp), allocatable :: jac_E(:,:,:,:,:)                                   ! jacobian of E(quilibrium) coordinate system
    real(dp), allocatable :: jac_F(:,:,:,:,:)                                   ! jacobian of F(lux) coordinate system with derivatives in the V(MEC) system
    real(dp), allocatable :: jac_FD(:,:,:,:,:)                                  ! jacobian of F(lux) coordinate system with derivatives in the F(lux) system
        
contains
    ! deallocates  metric  quantities  that  are  not  used  anymore  after  the
    ! equilibrium phase
    integer function dealloc_metric() result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'dealloc_metric'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate general variables
        deallocate(jac_F,h_F,g_F)
        deallocate(T_EF,T_FE)
        deallocate(det_T_EF,det_T_FE)
        deallocate(jac_E,g_E,h_E)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(T_VC,det_T_VC)
                deallocate(jac_C,g_C)
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    end function dealloc_metric
    
    ! deallocates  metric quantities  that are not  used anymore  after the
    ! calculation for a certain alpha
    subroutine dealloc_metric_final
        deallocate(g_FD,h_FD)
        deallocate(jac_FD)
    end subroutine dealloc_metric_final
end module metric_vars
