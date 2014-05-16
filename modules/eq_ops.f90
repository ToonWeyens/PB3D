!------------------------------------------------------------------------------!
!   Calculates the equilibrium quantities, making use of the metric_ops,       !
!   eq_vars, etc                                                               !
!------------------------------------------------------------------------------!
module eq_ops
    use num_vars, only: pi, dp
    use output_ops, only: print_ar_2, lvl_ud, writo
    
    implicit none
    private
    public calc_eq
    
    ! the possible derivatives of order i
    integer, allocatable :: derivs_0(:,:)                                       ! all possible derivatives of order 0
    integer, allocatable :: derivs_1(:,:)                                       ! all possible derivatives of order 1
    integer, allocatable :: derivs_2(:,:)                                       ! all possible derivatives of order 2
    integer, allocatable :: derivs_3(:,:)                                       ! all possible derivatives of order 3
    integer, allocatable :: derivs_4(:,:)                                       ! all possible derivatives of order 3

contains
    ! calculate the equilibrium quantities on a grid determined by straight field
    ! lines.
    subroutine calc_eq(alpha)
        use eq_vars, only: eqd_mesh, calc_mesh, calc_flux_q, &
            &check_mesh, init_eq, calc_RZL
        use B_vars, only: calc_B_V
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, &
            &calc_jac_V, calc_jac_F, calc_f_deriv, &
            &T_VF, T_FV, g_F, h_F, det_T_VF, det_T_FV, jac_F
        
        ! input / output
        real(dp) :: alpha
        
        call writo('Start setting up equilibrium quantities')
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call lvl_ud(1)
            call writo('Start initializing variables')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! initialize equilibrium quantities
            call writo('Initialize equilibrium quantities...')
            call init_eq
            
            ! initialize metric quantities
            call writo('Initialize metric quantities...')
            call init_metric
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Variables initialized')
            
            call writo('Start determining the equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! calculate mesh points (theta, zeta) that follow the magnetic field
            ! line
            call calc_mesh(alpha)
            
            ! check whether the mesh has been calculated correctly
            !call check_mesh(alpha)
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            !  calculate all  the  possible derivatives,  as  arguments for  the
            ! following subroutine
            call calc_derivs
            
            ! calculate  the   cylindrical  variables   R,  Z  and   lambda  and
            ! derivatives
            call writo('Calculate R,Z,L...')
            call calc_RZL(derivs_0)                                             ! order 0
            call calc_RZL(derivs_1)                                             ! order 1
            call calc_RZL(derivs_2)                                             ! order 2
            call calc_RZL(derivs_3)                                             ! order 3
            call calc_RZL(derivs_4)                                             ! order 4
            
            ! calculate flux quantities
            call calc_flux_q
            
            ! calculate the metrics in the cylindrical coordinate system
            call writo('Calculate g_C...')                                      ! h_C is not necessary
            call calc_g_C(derivs_0)                                             ! order 0
            call calc_g_C(derivs_1)                                             ! order 1
            call calc_g_C(derivs_2)                                             ! order 2
            call calc_g_C(derivs_3)                                             ! order 3
            
            ! calculate the jacobian in the cylindrical coordinate system
            call writo('Calculate jac_C...')
            call calc_jac_C(derivs_0)                                           ! order 0
            call calc_jac_C(derivs_1)                                           ! order 1
            call calc_jac_C(derivs_2)                                           ! order 2
            call calc_jac_C(derivs_3)                                           ! order 3
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call writo('Calculate T_VC...')
            call calc_T_VC(derivs_0)                                            ! order 0
            call calc_T_VC(derivs_1)                                            ! order 1
            call calc_T_VC(derivs_2)                                            ! order 2
            call calc_T_VC(derivs_3)                                            ! order 3
            
            ! calculate the metric factors in the VMEC coordinate system
            call writo('Calculate g_V...')
            call calc_g_V(derivs_0)                                             ! order 0
            call calc_g_V(derivs_1)                                             ! order 1
            call calc_g_V(derivs_2)                                             ! order 2
            call calc_g_V(derivs_3)                                             ! order 3
            
            ! calculate the jacobian in the VMEC coordinate system
            call writo('Calculate jac_V...')
            call calc_jac_V(derivs_0)                                           ! order 0
            call calc_jac_V(derivs_1)                                           ! order 1
            call calc_jac_V(derivs_2)                                           ! order 2
            call calc_jac_V(derivs_3)                                           ! order 3
            
            ! calculate the transformation matrix V(mec) -> F(lux)
            call writo('Calculate T_VF...')
            call calc_T_VF(derivs_0)                                            ! order 0
            call calc_T_VF(derivs_1)                                            ! order 1
            call calc_T_VF(derivs_2)                                            ! order 2
            call calc_T_VF(derivs_3)                                            ! order 3
            
            ! calculate the inverse of the transformation matrix T_VF
            call writo('Calculate T_FV...')
            call calc_inv_met(T_FV,T_VF,derivs_0)                               ! order 0
            call calc_inv_met(T_FV,T_VF,derivs_1)                               ! order 1
            call calc_inv_met(T_FV,T_VF,derivs_2)                               ! order 2
            call calc_inv_met(T_FV,T_VF,derivs_3)                               ! order 3
            call calc_inv_met(det_T_FV,det_T_VF,derivs_0)                       ! order 0
            call calc_inv_met(det_T_FV,det_T_VF,derivs_1)                       ! order 1
            call calc_inv_met(det_T_FV,det_T_VF,derivs_2)                       ! order 2
            call calc_inv_met(det_T_FV,det_T_VF,derivs_3)                       ! order 3
            
            ! calculate the metric factors in the Flux coordinate system
            call writo('Calculate g_F...')
            call calc_g_F(derivs_0)                                             ! order 0
            call calc_g_F(derivs_1)                                             ! order 1
            call calc_g_F(derivs_2)                                             ! order 2
            call calc_g_F(derivs_3)                                             ! order 3
            
            ! calculate the inverse of the metric factors g_F
            call writo('Calculate h_F...')
            call calc_inv_met(h_F,g_F,derivs_0)                                 ! order 0
            call calc_inv_met(h_F,g_F,derivs_1)                                 ! order 1
            call calc_inv_met(h_F,g_F,derivs_2)                                 ! order 2
            call calc_inv_met(h_F,g_F,derivs_3)                                 ! order 3
            
            ! calculate the jacobian in the Flux coordinate system
            call writo('Calculate jac_F...')
            call calc_jac_F(derivs_0)                                           ! order 0
            call calc_jac_F(derivs_1)                                           ! order 1
            call calc_jac_F(derivs_2)                                           ! order 2
            call calc_jac_F(derivs_3)                                           ! order 3
            
            ! calculate the derivatives in Flux coordinates from the derivatives
            ! in VMEC coordinates
            call writo('Calculate derivatives in flux coordinates...')
            call calc_f_deriv(g_F,derivs_1)                                     ! order 1
            call calc_f_deriv(g_F,derivs_2)                                     ! order 2
            call calc_f_deriv(g_F,derivs_3)                                     ! order 3
            call calc_f_deriv(h_F,derivs_1)                                     ! order 1
            call calc_f_deriv(h_F,derivs_2)                                     ! order 2
            call calc_f_deriv(h_F,derivs_3)                                     ! order 3
            call calc_f_deriv(jac_F,derivs_1)                                   ! order 1
            call calc_f_deriv(jac_F,derivs_2)                                   ! order 2
            call calc_f_deriv(jac_F,derivs_3)                                   ! order 3
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on equilibrium grid')
        
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end subroutine
    
    ! calculate all possible combinations of derivatives of a certain order
    subroutine calc_derivs
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        integer :: ci, cj, ck, cl                                               ! counters
        
        ! deallocate if needed
        if (allocated(derivs_0)) deallocate(derivs_0)
        if (allocated(derivs_1)) deallocate(derivs_1)
        if (allocated(derivs_2)) deallocate(derivs_2)
        if (allocated(derivs_3)) deallocate(derivs_3)
        if (allocated(derivs_4)) deallocate(derivs_4)
        allocate(derivs_0(3,1))
        allocate(derivs_1(3,3))
        allocate(derivs_2(3,6))
        allocate(derivs_3(3,10))
        allocate(derivs_4(3,15))
        
        ci = 1
        cj = 1
        ck = 1
        cl = 1
        
        derivs_0 = 0
        derivs_1 = 0
        derivs_2 = 0
        derivs_3 = 0
        derivs_4 = 0
        
        do id = 1,3
            derivs_1(id,ci) = derivs_1(id,ci) + 1
            ci = ci+1
            do jd = 1,id
                derivs_2(id,cj) = derivs_2(id,cj) + 1
                derivs_2(jd,cj) = derivs_2(jd,cj) + 1
                cj = cj+1
                do kd = 1,jd
                    derivs_3(id,ck) = derivs_3(id,ck) + 1
                    derivs_3(jd,ck) = derivs_3(jd,ck) + 1
                    derivs_3(kd,ck) = derivs_3(kd,ck) + 1
                    ck = ck+1
                    do ld = 1,kd
                        derivs_4(id,cl) = derivs_4(id,cl) + 1
                        derivs_4(jd,cl) = derivs_4(jd,cl) + 1
                        derivs_4(kd,cl) = derivs_4(kd,cl) + 1
                        derivs_4(ld,cl) = derivs_4(ld,cl) + 1
                        cl = cl+1
                    end do
                end do
            end do
        end do
    end subroutine
end module eq_ops
