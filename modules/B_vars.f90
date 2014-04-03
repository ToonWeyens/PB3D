!------------------------------------------------------------------------------!
!   This module contains operations concerning the magnetic field in various   !
!   coordinate systems                                                         !
!------------------------------------------------------------------------------!
module B_vars
    use num_vars, only: dp, pi
    use output_ops, only: print_ar_2
    
    implicit none
    private
    public calc_B_V, calc_B_F, &
        &B_V_sub, B_V_sub_H, B_F_sub, B_V, B_V_H, B_F, B_F_H
    
    real(dp), allocatable :: B_V_sub(:,:,:,:)                                   ! contravar. comp. of B in V(mec) coords, last index: (r,theta,zeta) (FM)
    real(dp), allocatable :: B_V_sub_H(:,:,:,:)                                 ! contravar. comp. of B in V(mec) coords, last index: (r,theta,zeta) (HM)
    real(dp), allocatable :: B_F_sub(:,:,:,:)                                   ! contravar. comp. of B in F(lux) coords, last index: (alpha,psi,zeta) (FM)
    real(dp), allocatable :: B_F_sub_H(:,:,:,:)                                 ! contravar. comp. of B in F(lux) coords, last index: (alpha,psi,zeta) (HM)
    real(dp), allocatable :: B_V(:,:), B_V_H(:,:)                               ! magnitude of B in V(MEC) coordinates (FM) and (HM)
    real(dp), allocatable :: B_F(:,:), B_F_H(:,:)                               ! magnitude of B in F(lux) coordinates (FM) and (HM)

contains
    ! Calculate  the components  of the  magnetic field  in the  VMEC coordinate
    ! system. Only the covariant components are calculated.
    subroutine calc_B_V
        use fourier_ops, only: mesh_cs, f2r
        use VMEC_vars, only: ntor, mpol, n_r, &
            &B_V_sub_c_M, B_V_sub_s_M, &                                        ! Coeff. of B_i in (co)sine series (last index: r,theta,phi) (FM, HM, HM)
            &B_V_c_H, B_V_s_H                                                   ! Coeff. of magnitude of B (HM)
        use eq_vars, only: n_par, theta, zeta, theta_H, zeta_H
        use utilities, only: calc_norm_deriv, f2h, h2f
        
        ! local variables
        real(dp) :: cs(0:mpol-1,-ntor:ntor,2)                                   ! (co)sines for all pol m and tor n
        integer :: id, jd, kd
        
        ! deallocate if allocated
        if (allocated(B_V_sub)) deallocate(B_V_sub)
        if (allocated(B_V_sub_H)) deallocate(B_V_sub_H)
        if (allocated(B_V)) deallocate(B_V)
        if (allocated(B_V_H)) deallocate(B_V_H)
        ! reallocate
        allocate(B_V_sub(n_par,n_r,4,3)); B_V_sub = 0.0_dp
        allocate(B_V_sub_H(n_par,n_r,4,3)); B_V_sub_H = 0.0_dp
        allocate(B_V(n_par,n_r)); B_V = 0.0_dp
        allocate(B_V_H(n_par,n_r)); B_V_H = 0.0_dp
        
        ! do calculations for all angular points and all components
        par: do id = 1, n_par                                                   ! parallel: along the magnetic field line
            ! calculate  the components  of the  variable B_V_sub  and their
            ! angular derivatives for all  normal points and current angular
            ! point
            perp: do kd = 1, n_r                                                ! perpendicular: normal to the flux surfaces
                ! FM quantities: r (s)
                cs = mesh_cs(mpol,ntor,theta(id,kd),zeta(id,kd))
                B_V_sub(id,kd,1,1) = f2r(B_V_sub_c_M(:,:,kd,1),&
                    &B_V_sub_s_M(:,:,kd,1),cs,mpol,ntor,1)
                B_V_sub(id,kd,3,1) = f2r(B_V_sub_c_M(:,:,kd,1),&
                    &B_V_sub_s_M(:,:,kd,1),cs,mpol,ntor,3)
                B_V_sub(id,kd,4,1) = f2r(B_V_sub_c_M(:,:,kd,1),&
                    &B_V_sub_s_M(:,:,kd,1),cs,mpol,ntor,4)
                
                ! HM quantities: theta (u), zeta (v)
                cs = mesh_cs(mpol,ntor,theta_H(id,kd),zeta_H(id,kd))
                comp: do jd = 2,3
                    B_V_sub_H(id,kd,1,jd) = f2r(B_V_sub_c_M(:,:,kd,jd),&
                        &B_V_sub_s_M(:,:,kd,jd),cs,mpol,ntor,1)
                    B_V_sub_H(id,kd,3,jd) = f2r(B_V_sub_c_M(:,:,kd,jd),&
                        &B_V_sub_s_M(:,:,kd,jd),cs,mpol,ntor,3)
                    B_V_sub_H(id,kd,4,jd) = f2r(B_V_sub_c_M(:,:,kd,jd),&
                        &B_V_sub_s_M(:,:,kd,jd),cs,mpol,ntor,4)
                end do comp
                B_V_H(id,kd) = f2r(B_V_c_H(:,:,kd),B_V_s_H(:,:,kd),cs,&
                    &mpol,ntor,1)
            end do perp
            
            ! numerically  calculate  normal  derivatives  at  the  currrent
            ! angular point
            B_V_sub(id,:,2,1) = calc_norm_deriv(B_V_sub(id,:,1,1),&
                &dfloat(n_r-1),.true.,.true.)
            
            comp3: do jd = 2,3                                                  ! components theta and phi (HM defined in VMEC)
                B_V_sub_H(id,:,2,jd) = calc_norm_deriv(B_V_sub_H(id,:,1,jd)&
                    &,dfloat(n_r-1),.false.,.false.)
            end do comp3
            
            ! convert from FM to HM: component r
            B_V_sub_H(id,:,1,1) = f2h(B_V_sub(id,:,1,1))
            B_V_sub_H(id,:,2,1) = calc_norm_deriv(B_V_sub(id,:,1,1),&
                    &dfloat(n_r-1),.true.,.false.)
            B_V_sub_H(id,:,3,1) = f2h(B_V_sub(id,:,1,1))
            B_V_sub_H(id,:,4,1) = f2h(B_V_sub(id,:,4,1))
            
            ! convert from HM to FM: components theta and zeta
            comp2: do jd = 2,3
                B_V_sub(id,:,1,jd) = h2f(B_V_sub_H(id,:,1,jd))
                B_V_sub(id,:,2,jd) = calc_norm_deriv(B_V_sub_H(id,:,1,jd),&
                        &n_r-1._dp,.false.,.true.)
                B_V_sub(id,:,3,jd) = h2f(B_V_sub_H(id,:,3,jd))
                B_V_sub(id,:,4,jd) = h2f(B_V_sub_H(id,:,4,jd))
            end do comp2
            ! and magnitude
            B_V(id,:) = h2f(B_V_H(id,:))
        end do par
    end subroutine
    
    ! Calculate  the components  of the  magnetic field  in the  flux coordinate
    ! system by transforming the VMEC components. In the flux coordinate systen,
    ! only  the  covariant  components  are  calculated,  as  the  contravariant
    ! components are trivial (e.g. B^alpha = B^phi = 0 and B^theta = 1/J_F)
    ! These are then compared to the calculated  components from B_i = e_i * B =
    ! e_theta*e_i/J_F = g_F(3,i)/J_F.
    subroutine calc_B_F
        use VMEC_vars, only: n_r
        use utilities, only: matvec_mult
        use metric_ops, only: V2F_dn, V2F_dn_H, jac_F, jac_F_H, g_F, g_F_H
        use eq_vars, only: n_par
        
        ! local variables
        integer :: id, kd
        real(dp) :: B_F_sub_H_alt(n_par,n_r,3)
        real(dp) :: B_F_sub_alt(n_par,n_r,3)
        real(dp) :: B_F_alt(n_par,n_r)
        real(dp) :: B_F_H_alt(n_par,n_r)
        
        ! deallocate if allocated
        if (allocated(B_F_sub)) deallocate(B_F_sub)
        if (allocated(B_F_sub_H)) deallocate(B_F_sub_H)
        if (allocated(B_F)) deallocate(B_F)
        if (allocated(B_F_H)) deallocate(B_F_H)
        ! reallocate
        allocate(B_F_sub(n_par,n_r,4,3)); B_F_sub = 0.0_dp
        allocate(B_F_sub_H(n_par,n_r,4,3)); B_F_sub_H = 0.0_dp
        allocate(B_F(n_par,n_r)); B_F = 0.0_dp
        allocate(B_F_H(n_par,n_r)); B_F_H = 0.0_dp
        
        ! transform from VMEC to flux components using the transformation matrix
        ! V2F_dn.  Be  careful with  the  derivatives:  Also take  into  account
        ! the  derivatives fo  the  unit vectors.  Also, make  sure  to use  the
        ! correct variant  (FM or  HM) for  the normal  derivatives, as  for the
        ! s-component, the FM  version is more accurate while for  the theta and
        ! zeta-components, the HM version is more accurate
        par: do id = 1, n_par                                                   ! parallel: along the magnetic field line
            ! calculate  the  values  of  the contravar.  components  of  B,  no
            ! derivatives
            perp: do kd = 1, n_r                                                ! perpendicular: normal to the flux surfaces
                B_F_sub(id,kd,1,:) = &
                    &matvec_mult(V2F_dn(:,:,id,kd),B_V_sub(id,kd,1,:))          ! B_F_sub = V2F_dn * B_V_sub (FM)
                B_F_sub_H(id,kd,1,:) = &
                    &matvec_mult(V2F_dn_H(:,:,id,kd),B_V_sub_H(id,kd,1,:))      ! B_F_sub = V2F_dn * B_V_sub (HM)
            end do perp
        end do par
        ! the magnitude of a vector stays constant:
        B_F = B_V
        B_F_H = B_V_H
        
        ! check this against an alternative calculation using g_F(3,i)/J_F
        perp2: do kd = 1, n_r
            par2: do id = 1, n_par
                B_F_sub_alt(id,kd,1) = g_F(3,1,id,kd) / jac_F(id,kd)
                B_F_sub_H_alt(id,kd,2) = g_F_H(3,2,id,kd) / jac_F_H(id,kd)
                B_F_sub_H_alt(id,kd,3) = g_F_H(3,3,id,kd) / jac_F_H(id,kd)
            end do par2
        end do perp2
        
        write(*,*) 'B_F_sub (alpha) = '
        call print_ar_2(B_F_sub(:,:,1,1))
        write(*,*) 'B_F_sub_alt (alpha) = '
        call print_ar_2(B_F_sub_alt(:,:,1))
        write(*,*) 'B_F_sub_H (psi) = '
        call print_ar_2(B_F_sub_H(:,:,1,2))
        write(*,*) 'B_F_sub_alt_H (psi) = '
        call print_ar_2(B_F_sub_H_alt(:,:,2))
        write(*,*) 'B_F_sub_H (theta) = '
        call print_ar_2(B_F_sub_H(:,:,1,3))
        write(*,*) 'B_F_sub_alt_H (theta) = '
        call print_ar_2(B_F_sub_H_alt(:,:,3))
        write(*,*) ''
        read(*,*)
        
        do kd = 1, n_r
            do id = 1, n_par
                B_F_alt(id,kd) = sqrt(B_F_sub(id,kd,1,3) / jac_F(id,kd) )
                B_F_H_alt(id,kd) = sqrt(B_F_sub_H(id,kd,1,3) / jac_F_H(id,kd) )
            end do
        end do
        
        write(*,*) 'magn of B = '
        call print_ar_2(B_V)
        write(*,*) 'magn of B, alternative = '
        call print_ar_2(B_F_alt(:,:))
        write(*,*) 'magn of B_H = '
        call print_ar_2(B_V_H)
        write(*,*) 'magn of B, alternative = '
        call print_ar_2(B_F_H_alt(:,:))
        read(*,*)
        
        !read(*,*)
        !write(*,*) 'B_F_sub at, alpha component, at kd = 5 :'
        !call print_ar_2(B_F_sub(:,5,:,1))
        !write(*,*) 'B_F_sub at, psi component, at kd = 5 :'
        !call print_ar_2(B_F_sub(:,5,:,2))
        !write(*,*) 'B_F_sub at, theta component, at kd = 5 :'
        !call print_ar_2(B_F_sub(:,5,:,3))
        !read(*,*)
        
        !do kd = 1, n_r
            !write(*,*) ' B_F_sub at kd =', kd
            !call print_ar_2(B_F_sub(:,kd,1,:))
            !write(*,*) ' B_F_sub_alt at kd =', kd
            !call print_ar_2(B_F_sub_alt(:,kd,:))
            !read(*,*)
        !end do
        
        !! alternative calculation 2 
        !perp3: do kd = 1, n_r
            !par3: do id = 1, n_par
                !B_V_sub_alt(:,:,1) = flux_p(kd,2)/(2*pi) * 1 / jac_V(id,kd) * &
                    !&((1-q_saf(kd,1)*lam_H(id,kd,4))*R(id,kd,3) + &
                    !&q_saf(kd,1)*(1+lam_H(id,kd,3))*R(id,kd,4))
                !B_V_sub_alt(:,:,2) = flux_p(kd,2)/(2*pi) * 1 / jac_V(id,kd) * &
                    !&(-q_saf(kd,1) * (1 - lam_H(id,kd,3)))
                !B_V_sub_alt(:,:,3) = flux_p(kd,2)/(2*pi) * 1 / jac_V(id,kd) * &
                    !&((1-q_saf(kd,1)*lam_H(id,kd,4))*Z(id,kd,3) + &
                    !&q_saf(kd,1)*(1+lam_H(id,kd,3))*Z(id,kd,4))
            !end do par3
        !end do perp3
        !write(*,*) 'THE VARIABLE (NO DERIVATIVES)'
        !write(*,*) ' B_V_sub: r component'
        !call print_ar_2(B_V_sub(:,:,1,1))
        !write(*,*) 'B_V_sub_alt, r component: '
        !call print_ar_2(B_V_sub_alt(:,:,1))
        !write(*,*) ' B_V_sub: theta component'
        !call print_ar_2(B_V_sub(:,:,1,2))
        !write(*,*) 'B_V_sub_alt, theta component: '
        !call print_ar_2(B_V_sub_alt(:,:,2))
        !write(*,*) ' B_V_sub: zeta component'
        !call print_ar_2(B_V_sub(:,:,1,3))
        !write(*,*) 'B_V_sub_alt, zeta component: '
        !call print_ar_2(B_V_sub_alt(:,:,3))
        
        !write(*,*) 'R = '
        !call print_ar_2(R(:,:,1))
    end subroutine
end module B_vars
