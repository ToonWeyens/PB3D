!------------------------------------------------------------------------------!
!   Variables that have to do with the metric elements                         !
!------------------------------------------------------------------------------!
module metric_vars
#include <PB3D_macros.h>
    use num_vars, only: dp, max_deriv, max_str_ln
    use str_ops, only: r2str, i2str
    use output_ops, only: print_GP_3D
    
    implicit none
    private
    public create_metric, dealloc_metric, dealloc_metric_final, &
        &metric_type
    
    ! metric type
    ! The arrays here are of the form:
    !   - (angle_1,angle_2,r,Dr,Dangle_1,Dangle_2)      for normal quantities
    !   - (angle_1,angle_2,r,6/9,Dr,Dangle_1,Dangle_2)  for metric matrices
    ! where  a discussion  of the  coordinates is  given in  the description  of
    ! eq_type. Also, like  with the equilibrium type, the metric  type should be
    ! complemented by grid  type. The last index for  metric matrices correspond
    ! to the 9 or 6 (symmetric) different values:
    !   (1 4 7)      (1    )
    !   (2 5 8)  or  (2 4  )
    !   (3 6 9)      (3 5 6)
    ! Note that Fortran only allows for 7 dimensions in arrays.
    type :: metric_type
        ! upper (h) and lower (g) metric factors
        real(dp), allocatable :: g_C(:,:,:,:,:,:,:)                             ! in the C(ylindrical) coordinate system
        real(dp), allocatable :: g_E(:,:,:,:,:,:,:)                             ! in the E(quilibrium) coordinate system
        real(dp), allocatable :: h_E(:,:,:,:,:,:,:)                             ! in the E(quilibrium) coordinate system
        real(dp), allocatable :: g_F(:,:,:,:,:,:,:)                             ! in the F(lux) coordinate system with derivatves in the V(MEC) system
        real(dp), allocatable :: h_F(:,:,:,:,:,:,:)                             ! in the F(lux) coordinate system with derivatves in the V(MEC) system
        real(dp), pointer :: g_FD(:,:,:,:,:,:,:)                                ! in the F(lux) coordinate system with derivatives in the F(lux) system
        real(dp), pointer :: h_FD(:,:,:,:,:,:,:)                                ! in the F(lux) coordinate system with derivatives in the F(lux) system
        ! upper and lower transformation matrices
        real(dp), allocatable :: T_VC(:,:,:,:,:,:,:)                            ! C(ylindrical) to V(MEC) (lower)
        real(dp), allocatable :: T_EF(:,:,:,:,:,:,:)                            ! E(quilibrium) to F(lux) (upper)
        real(dp), allocatable :: T_FE(:,:,:,:,:,:,:)                            ! F(lux) to E(quilibrium) (lower)
        ! determinants of transformation matrices
        real(dp), allocatable :: det_T_VC(:,:,:,:,:,:)                          ! determinant of T_VC
        real(dp), allocatable :: det_T_EF(:,:,:,:,:,:)                          ! determinant of T_EF
        real(dp), allocatable :: det_T_FE(:,:,:,:,:,:)                          ! determinant of T_FE
        ! Jacobians
        real(dp), allocatable :: jac_C(:,:,:,:,:,:)                             ! jacobian of C(ylindrical) coordinate system
        real(dp), allocatable :: jac_E(:,:,:,:,:,:)                             ! jacobian of E(quilibrium) coordinate system
        real(dp), allocatable :: jac_F(:,:,:,:,:,:)                             ! jacobian of F(lux) coordinate system with derivatives in the V(MEC) system
        real(dp), pointer :: jac_FD(:,:,:,:,:,:)                                ! jacobian of F(lux) coordinate system with derivatives in the F(lux) system
    end type
    
contains
    ! initialize metric variables
    integer function create_metric(grid,met) result(ierr)
        use num_vars, only: eq_style
        use grid_vars, only: grid_type
        
        character(*), parameter :: rout_name = 'create_metric'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        type(grid_type), intent(in) :: grid                                     ! equilibrium grid
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: dims(3)                                                      ! dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! set up local variables
        dims = [grid%n(1),grid%n(2),grid%grp_n_r]
        
        ! initialize variables that are used for all equilibrium styles
        ! g_E
        allocate(met%g_E(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! h_E
        allocate(met%h_E(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! g_F
        allocate(met%g_F(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! h_F
        allocate(met%h_F(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! g_FD
        allocate(met%g_FD(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! h_FD
        allocate(met%h_FD(dims(1),dims(2),dims(3),6,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! jac_E
        allocate(met%jac_E(dims(1),dims(2),dims(3),&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! jac_F
        allocate(met%jac_F(dims(1),dims(2),dims(3),&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! jac_FD
        allocate(met%jac_FD(dims(1),dims(2),dims(3),&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! T_EF
        allocate(met%T_EF(dims(1),dims(2),dims(3),9,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! T_FE
        allocate(met%T_FE(dims(1),dims(2),dims(3),9,&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! det_T_EF
        allocate(met%det_T_EF(dims(1),dims(2),dims(3),&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! det_T_FE
        allocate(met%det_T_FE(dims(1),dims(2),dims(3),&
            &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
        
        ! initialize variables that are  specifici to which equilibrium style is
        ! being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                ! g_C
                allocate(met%g_C(dims(1),dims(2),dims(3),6,&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                
                ! T_VC
                allocate(met%T_VC(dims(1),dims(2),dims(3),9,&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                
                ! det_T_VC
                allocate(met%det_T_VC(dims(1),dims(2),dims(3),&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
                
                ! jac_C
                allocate(met%jac_C(dims(1),dims(2),dims(3),&
                    &0:max_deriv+1,0:max_deriv+1,0:max_deriv+1))
            case (2)                                                            ! HELENA
                ! nothing
            case default
                err_msg = 'No equilibrium style associated with '//&
                    &trim(i2str(eq_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
    contains
    end function create_metric
    
    ! deallocates  metric  quantities  that  are  not  used  anymore  after  the
    ! equilibrium phase
    integer function dealloc_metric(met) result(ierr)
        use num_vars, only: eq_style
        
        character(*), parameter :: rout_name = 'dealloc_metric'
        
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! deallocate general variables
        deallocate(met%h_E)
        deallocate(met%g_E)
        deallocate(met%h_F)
        deallocate(met%g_F)
        deallocate(met%T_EF)
        deallocate(met%T_FE)
        deallocate(met%det_T_EF,met%det_T_FE)
        deallocate(met%jac_F,met%jac_E)
        
        ! choose which equilibrium style is being used:
        !   1:  VMEC
        !   2:  HELENA
        select case (eq_style)
            case (1)                                                            ! VMEC
                deallocate(met%g_C)
                deallocate(met%T_VC)
                deallocate(met%det_T_VC)
                deallocate(met%jac_C)
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
    subroutine dealloc_metric_final(met)
        ! input / output
        type(metric_type), intent(inout) :: met                                 ! metric to be created
        
        ! deallocate variables
        nullify(met%g_FD)
        nullify(met%h_FD)
        nullify(met%jac_FD)
    end subroutine dealloc_metric_final
end module metric_vars
