!------------------------------------------------------------------------------!
!   Calculates the equilibrium quantities, making use of the metric_ops,       !
!   eq_vars, etc                                                               !
!------------------------------------------------------------------------------!
module eq_ops
#include <PB3D_macros.h>
    use num_vars, only: pi, dp
    use output_ops, only: print_ar_2, lvl_ud, writo
    
    implicit none
    private
    public calc_eq

contains
    ! calculate the equilibrium quantities on a grid determined by straight field
    ! lines.
    integer function calc_eq(alpha) result(ierr)
        use eq_vars, only: calc_mesh, calc_flux_q, dealloc_eq_vars, &
            &check_mesh, init_eq, calc_RZL, q_saf, q_saf_FD, flux_p, flux_p_FD,&
            &flux_t, flux_t_FD, pres, pres_FD
        use metric_ops, only: calc_g_C, calc_g_C, calc_T_VC, calc_g_V, &
            &init_metric, calc_T_VF, calc_inv_met, calc_g_F, calc_jac_C, &
            &calc_jac_V, calc_jac_F, calc_f_deriv, dealloc_metric_vars, &
            &T_VF, T_FV, g_F, h_F, det_T_VF, det_T_FV, jac_F, g_F_FD, h_F_FD, &
            &jac_F_FD
        use utilities, only: derivs
        use num_vars, only: max_deriv
        
        character(*), parameter :: rout_name = 'calc_eq'
        
        ! input / output
        real(dp) :: alpha
        
        ! local variables
        integer :: id
        
        ! initialize ierr
        ierr = 0
        
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
            ierr = calc_mesh(alpha)
            CHCKERR('')
            
            ! check whether the mesh has been calculated correctly
            ! if so, limit the normal range of theta and zeta to the local range
            ierr = check_mesh(alpha)
            CHCKERR('')
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium grid determined')
            
            call writo('Calculating equilibrium quantities on equilibrium grid')
            ! 2----------------------------------------------------------------
            call lvl_ud(1)
            
            ! calculate  the   cylindrical  variables   R,  Z  and   lambda  and
            ! derivatives
            call writo('Calculate R,Z,L...')
            do id = 0,2
                ierr = calc_RZL(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate flux quantities
            ierr = calc_flux_q()
            CHCKERR('')
            
            ! calculate the metrics in the cylindrical coordinate system
            call writo('Calculate g_C...')                                      ! h_C is not necessary
            do id = 0,1
                ierr = calc_g_C(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the jacobian in the cylindrical coordinate system
            call writo('Calculate jac_C...')
            do id = 0,1
                ierr = calc_jac_C(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the transformation matrix C(ylindrical) -> V(mec)
            call writo('Calculate T_VC...')
            do id = 0,1
                ierr = calc_T_VC(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the metric factors in the VMEC coordinate system
            call writo('Calculate g_V...')
            do id = 0,1
                ierr = calc_g_V(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the jacobian in the VMEC coordinate system
            call writo('Calculate jac_V...')
            do id = 0,1
                ierr = calc_jac_V(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the transformation matrix V(mec) -> F(lux)
            call writo('Calculate T_VF...')
            do id = 0,1
                ierr = calc_T_VF(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the inverse of the transformation matrix T_VF
            call writo('Calculate T_FV...')
            do id = 0,1
                ierr = calc_inv_met(T_FV,T_VF,derivs(id))
                CHCKERR('')
                ierr = calc_inv_met(det_T_FV,det_T_VF,derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the metric factors in the Flux coordinate system
            call writo('Calculate g_F...')
            do id = 0,1
                ierr = calc_g_F(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the inverse h_F of the metric factors g_F
            call writo('Calculate h_F...')
            do id = 0,1
                ierr = calc_inv_met(h_F,g_F,derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the jacobian in the Flux coordinate system
            call writo('Calculate jac_F...')
            do id = 0,1
                ierr = calc_jac_F(derivs(id))
                CHCKERR('')
            end do
            
            ! calculate the derivatives in Flux coordinates from the derivatives
            ! in VMEC coordinates
            call writo('Calculate derivatives in flux coordinates...')
            do id = 0,1
                ierr = calc_f_deriv(g_F,T_FV,g_F_FD,max_deriv-[1,1,1],&
                    &derivs(id))                                                ! g_F
                CHCKERR('')
                ierr = calc_f_deriv(h_F,T_FV,h_F_FD,max_deriv-[1,1,1],&
                    &derivs(id))                                                ! h_F
                CHCKERR('')
                ierr = calc_f_deriv(jac_F,T_FV,jac_F_FD,max_deriv-[1,1,1],&
                    &derivs(id))                                                ! jac_F
                CHCKERR('')
                ierr = calc_f_deriv(q_saf,T_FV(1,:,2,1,:,0,0),q_saf_FD(:,id),&
                    &max_deriv(1)-1,id)                                         ! q_saf
                CHCKERR('')
                ierr = calc_f_deriv(flux_p,T_FV(1,:,2,1,:,0,0),flux_p_FD(:,id),&
                    &max_deriv(1)-1,id)                                         ! flux_p
                CHCKERR('')
                ierr = calc_f_deriv(flux_t,T_FV(1,:,2,1,:,0,0),flux_t_FD(:,id),&
                    &max_deriv(1)-1,id)                                         ! flux_t
                CHCKERR('')
                ierr = calc_f_deriv(pres,T_FV(1,:,2,1,:,0,0),pres_FD(:,id),&
                    &max_deriv(1)-1,id)                                         ! pres
                CHCKERR('')
            end do
            
            ! deallocate unusd equilibrium quantities
            call writo('Deallocate unused equilibrium and metric quantities...')
            call dealloc_eq_vars
            call dealloc_metric_vars
            
            call lvl_ud(-1)
            ! 2----------------------------------------------------------------
            call writo('Equilibrium quantities calculated on equilibrium grid')
        
        call lvl_ud(-1)
        ! 1--------------------------------------------------------------------
        ! 1--------------------------------------------------------------------
        call writo('Done setting up equilibrium quantities')
    end function calc_eq
end module eq_ops
