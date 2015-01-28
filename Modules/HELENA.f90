!------------------------------------------------------------------------------!
!   operations and variable pertaining to the HELENA equilibrim                !
!------------------------------------------------------------------------------!
module HELENA
#include <PB3D_macros.h>
    use messages, only: writo, lvl_ud
    use num_vars, only: dp, max_str_ln
    use output_ops, only: print_GP_2D, print_GP_3D
    use str_ops, only: i2str, r2str
    
    implicit none
    private
    public read_HEL, dealloc_HEL, dealloc_HEL_final, &
        &R_0_H, B_0_H, p0, qs, flux_H, nchi, chi_H, ias, h_H_11, h_H_12, &
        &h_H_33, RBphi, R_H, Z_H, h_H_11_full, h_H_12_full, h_H_33_full
    
    ! global variables
    real(dp) :: R_0_H = 1.0_dp                                                  ! R of magnetic axis (normalization constant)
    real(dp) :: B_0_H = 1.0_dp                                                  ! B at magnetic axis (normalization constant)
    real(dp), allocatable :: chi_H(:)                                           ! poloidal angle
    real(dp), allocatable :: flux_H(:)                                          ! normal coordinate values
    real(dp), allocatable :: p0(:)                                              ! pressure profile
    real(dp), allocatable :: qs(:)                                              ! safety factor
    real(dp), allocatable :: RBphi(:)                                           ! R B_phi
    real(dp), allocatable :: h_H_11_full(:,:)                                   ! full upper metric factor 11 (gem11)
    real(dp), allocatable :: h_H_12_full(:,:)                                   ! full upper metric factor 12 (gem12)
    real(dp), allocatable :: h_H_33_full(:,:)                                   ! full upper metric factor 33 (1/gem33)
    real(dp), allocatable :: h_H_11(:,:)                                        ! adapted upper metric factor 11 (gem11) (see adapt_HEL_to_eq)
    real(dp), allocatable :: h_H_12(:,:)                                        ! adapted upper metric factor 12 (gem12) (see adapt_HEL_to_eq)
    real(dp), allocatable :: h_H_33(:,:)                                        ! adapted upper metric factor 33 (1/gem33) (see adapt_HEL_to_eq)
    real(dp), allocatable :: R_H(:,:)                                           ! major radius R (xout)
    real(dp), allocatable :: Z_H(:,:)                                           ! height Z (yout)
    integer :: nchi                                                             ! nr. of poloidal points (nchi)
    integer :: ias                                                              ! 0 if top-bottom symmetric, 1 if not

contains
    ! Reads the HELENA equilibrium data
    ! (from HELENA routine IODSK)
    ! Note: The variables in HELENA output are normalized wrt.:
    !   - R_m: radius of magnetic axis at equilibrium
    !   - B_m: magnetic field at magnetic axis at equilibrium ,
    ! These have to be specified (see global variables above)
    ! However, the variables R_H and Z_H are also normalized wrt.:
    !   - R_H = (R-R_0)/a
    !   - Z_H = Z/a ,
    ! where R_0 and a are found through the variable "radius" and "eps":
    !   - radius = a / R_m
    !   - eps = a / R_0 ,
    ! which link R_0 and a to R_m.
    ! Furthermore, the varaible "cs" contains the sqrt of the normalized flux on
    ! the  normal positions,  where  the normalization  factor  is contained  in
    ! "cpsurf", which is the poloidal flux at the surface.
    ! Finally, some variables are not tabulated on the magnetic axis.
    ! [MPI] only global master
    integer function read_HEL(n_r_eq,use_pol_flux_H) result(ierr)
        use num_vars, only: glb_rank, eq_name, eq_i, mu_0
        
        character(*), parameter :: rout_name = 'read_HEL'
        
        ! input / output
        integer, intent(inout) :: n_r_eq                                        ! nr. of normal points in equilibrium grid
        logical, intent(inout) :: use_pol_flux_H                                ! .true. if HELENA equilibrium is based on pol. flux
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: dqs(:)                                         ! derivative of q
        real(dp), allocatable :: curj(:)                                        ! toroidal current
        real(dp), allocatable :: vx(:), vy(:)                                   ! R and Z of surface
        real(dp) :: dj0, dje                                                    ! derivative of toroidal current on axis and surface
        real(dp) :: dp0, dpe                                                    ! normal derivative of pressure on axis and surface
        real(dp) :: dRBphi0, dRBphie                                            ! normal derivative of R B_phi on axis and surface
        real(dp) :: radius, raxis                                               ! minor radius, major radius
        real(dp) :: eps                                                         ! inverse aspect ratio
        real(dp) :: cpsurf                                                      ! poloidal flux on surface
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then                                                 ! only global master
            call writo('Reading data from HELENA output "' &
                &// trim(eq_name) // '"')
            call lvl_ud(1)
            
            ! set error messages
            err_msg = 'Failed to read HELENA output file'
            
            ! rewind input file
            rewind(eq_i)
            
            ! read mapped equilibrium from disk
            read(eq_i,*,IOSTAT=ierr) n_r_eq,ias                                 ! nr. normal points (JS0), top-bottom symmetry
            CHCKERR(err_msg)
            n_r_eq = n_r_eq + 1                                                 ! HELENA outputs nr. of normal points - 1
            
            allocate(flux_H(n_r_eq))                                            ! flux, derived from normal coordinate
            read(eq_i,*,IOSTAT=ierr) (flux_H(kd),kd=1,n_r_eq)                   ! it is squared below, after reading cpsurf
            CHCKERR(err_msg)
            
            allocate(qs(n_r_eq))                                                ! safety factor
            read(eq_i,*,IOSTAT=ierr) (qs(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            
            allocate(dqs(n_r_eq))                                               ! derivative of safety factor
            read(eq_i,*,IOSTAT=ierr) dqs(1),dqs(n_r_eq)                         ! first point, last point
            CHCKERR(err_msg)
            
            read(eq_i,*,IOSTAT=ierr) (dqs(kd),kd=2,n_r_eq)                      ! second to last point (again)
            CHCKERR(err_msg)
            
            allocate(curj(n_r_eq))                                              ! toroidal current
            read(eq_i,*,IOSTAT=ierr) (curj(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            
            read(eq_i,*,IOSTAT=ierr) dj0,dje                                    ! derivative of toroidal current at first and last point
            CHCKERR(err_msg)
            
            read(eq_i,*,IOSTAT=ierr) nchi                                       ! nr. poloidal points
            CHCKERR(err_msg)
            
            allocate(chi_H(nchi))                                               ! poloidal points
            read(eq_i,*,IOSTAT=ierr) (chi_H(id),id=1,nchi)
            CHCKERR(err_msg)
            
            allocate(h_H_11_full(nchi,n_r_eq))                                  ! upper metric factor 1,1
            read(eq_i,*,IOSTAT=ierr) &
                &(h_H_11_full(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,&
                &n_r_eq*nchi)                                                   ! (gem11)
            CHCKERR(err_msg)
            h_H_11_full(:,1) = 0._dp                                            ! first normal point is not given, so set to zero
            h_H_11_full = h_H_11_full * (R_0_H * B_0_H)**2                      ! rescale h_H_11_full
            
            allocate(h_H_12_full(nchi,n_r_eq))                                  ! upper metric factor 1,2
            read(eq_i,*,IOSTAT=ierr) &
                &(h_H_12_full(mod(id-1,nchi)+1,(id-1)/nchi+1),&
                &id=nchi+1,n_r_eq*nchi)                                         ! (gem12)
            CHCKERR(err_msg)
            h_H_12_full(:,1) = 0._dp                                            ! first normal point is not given, so set to zero
            h_H_12_full = h_H_12_full * B_0_H                                   ! rescale h_H_12_full
            
            read(eq_i,*,IOSTAT=ierr) cpsurf, radius                             ! poloidal flux on surface, minor radius
            CHCKERR(err_msg)
            cpsurf = cpsurf * R_0_H**2 * B_0_H                                  ! back to real space
            flux_H = flux_H**2 * cpsurf                                         ! rescale poloidal flux
            
            allocate(h_H_33_full(nchi,n_r_eq))                                  ! upper metric factor 3,3
            read(eq_i,*,IOSTAT=ierr) &
                &(h_H_33_full(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,&
                &n_r_eq*nchi)                                                   ! (gem33)
            h_H_33_full(:,:) = 1._dp/h_H_33_full(:,:)                           ! HELENA gives R^2, but need 1/R^2
            h_H_33_full(:,1) = 0._dp                                            ! first normal point is not given, so set to zero
            h_H_11_full = h_H_11_full / (R_0_H**2)                              ! rescale h_H_33_full
            
            read(eq_i,*,IOSTAT=ierr) raxis                                      ! major radius
            CHCKERR(err_msg)
            
            allocate(p0(n_r_eq))                                                ! pressure profile
            read(eq_i,*,IOSTAT=ierr) (p0(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            p0 = p0 * B_0_H**2/mu_0                                             ! rescale pressure
            
            read(eq_i,*,IOSTAT=ierr) dp0,dpe                                    ! derivarives of pressure on axis and surface
            CHCKERR(err_msg)
            
            allocate(RBphi(n_r_eq))                                             ! R B_phi (= F)
            read(eq_i,*,IOSTAT=ierr) (RBphi(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            RBphi = RBphi * R_0_H * B_0_H                                       ! rescale F
            
            read(eq_i,*,IOSTAT=ierr) dRBphi0,dRBphie                            ! derivatives of R B_phi on axis and surface
            CHCKERR(err_msg)
            
            allocate(vx(nchi))                                                  ! R B_phi
            read(eq_i,*,IOSTAT=ierr) (vx(id),id=1,nchi)                         ! R on surface
            CHCKERR(err_msg)
            
            allocate(vy(nchi))                                                  ! R B_phi
            read(eq_i,*,IOSTAT=ierr) (vy(id),id=1,nchi)                         ! Z on surface
            CHCKERR(err_msg)
            
            read(eq_i,*,IOSTAT=ierr) eps                                        ! inerse aspect ratio
            CHCKERR(err_msg)
            
            allocate(R_H(nchi,n_r_eq))                                          ! major radius R
            R_H(:,1) = 0._dp                                                    ! values on axis are not given by HELENA -> set to zero
            read(eq_i,*,IOSTAT=ierr) &
                &(R_H(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)    ! (xout)
            CHCKERR(err_msg)
            R_H = R_0_H*radius*(1._dp/eps + R_H)                                ! back to real space
            
            allocate(Z_H(nchi,n_r_eq))                                          ! height Z
            Z_H(:,1) = 0._dp                                                    ! values on axis are not given by HELENA -> set to zero
            read(eq_i,*,IOSTAT=ierr) &
                &(Z_H(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)    ! (yout)
            CHCKERR(err_msg)
            Z_H = R_0_H*radius*Z_H                                              ! back to real space
           
            ! HELENA always uses the poloidal flux
            use_pol_flux_H = .true.
            
            call writo('HELENA output given on '//trim(i2str(nchi))//&
                &' poloidal and '//trim(i2str(n_r_eq))//' normal points')
            call lvl_ud(-1)
            call writo('Grid parameters successfully read')
            
            ! test whether the metric quantities are correctly derived
            call writo('Checking consinstency of metric factors')
            
            call lvl_ud(1)
            ierr = check_metrics(n_r_eq)
            CHCKERR('')
            call lvl_ud(-1)
        end if
    end function read_HEL
    
    ! checks whether the metric elements  provided by HELENA are consistent with
    ! a direct calculation using the coordinate transformations:
    !   |nabla psi|^2       = 1/jac^2 ((dZ/dchi)^2 + (dR/dchi)^2)
    !   nabla psi nabla chi = 1/jac^2 (dZ/dchi dZ/dpsi + dR/dchi dR/dpsi)
    !   |nabla chi|^2       = 1/jac^2 ((dZ/dpsi)^2 + (dR/dpsi)^2)
    !   |nabla phi|^2       = 1/R^2
    ! with jac = dZ/dpsi dR/dchi - dR/dpsi dZ/dchi
    integer function check_metrics(n_r) result(ierr)
        use utilities, only: calc_deriv
        
        character(*), parameter :: rout_name = 'check_metrics'
        
        ! input / output
        integer, intent(inout) :: n_r                                           ! nr. of normal points in equilibrium grid
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: Rchi(:,:), Rpsi(:,:)                           ! chi and psi derivatives of R
        real(dp), allocatable :: Zchi(:,:), Zpsi(:,:)                           ! chi and psi derivatives of Z
        real(dp), allocatable :: jac(:,:)                                       ! jac as defined above
        real(dp), allocatable :: h_H_11_alt(:,:)                                ! alternative calculation for upper metric factor 11
        real(dp), allocatable :: h_H_12_alt(:,:)                                ! alternative calculation for upper metric factor 12
        real(dp), allocatable :: h_H_33_alt(:,:)                                ! alternative calculation for upper metric factor 33
        
        ! initialize ierr
        ierr = 0
        
        ! calculate  the   auxiliary  quantities  Zchi,  zpsi,   Rchi  and  Rpsi
        ! containing the derivatives as well as jac
        allocate(Rchi(nchi,n_r),Rpsi(nchi,n_r))
        allocate(Zchi(nchi,n_r),Zpsi(nchi,n_r))
        allocate(jac(nchi,n_r))
        
        do id = 1,nchi
            ierr = calc_deriv(R_H(id,:),Rpsi(id,:),flux_H,1,1)
            CHCKERR('')
            ierr = calc_deriv(Z_H(id,:),Zpsi(id,:),flux_H,1,1)
            CHCKERR('')
        end do
        do kd = 1,n_r
            ierr = calc_deriv(R_H(:,kd),Rchi(:,kd),chi_H,1,1)
            CHCKERR('')
            ierr = calc_deriv(Z_H(:,kd),Zchi(:,kd),chi_H,1,1)
            CHCKERR('')
        end do
        jac = Zpsi*Rchi - Zchi*Rpsi
        
        ! calculate the metric factors directly
        allocate(h_H_11_alt(nchi,n_r))
        allocate(h_H_12_alt(nchi,n_r))
        allocate(h_H_33_alt(nchi,n_r))
        
        h_H_11_alt = 1._dp/(jac**2) * (Zchi**2 + Rchi**2)
        h_H_12_alt = -1._dp/(jac**2) * (Zchi*Zpsi + Rchi*Rpsi)
        h_H_33_alt = 1._dp/(R_H**2)
        
        ! output results
        call writo('maximum error in h_H_11 = '//&
            &trim(r2str(maxval(abs(h_H_11_full(:,3:n_r)-h_H_11_alt(:,3:n_r))))))
        call writo('maximum error in h_H_12 = '//&
            &trim(r2str(maxval(abs(h_H_12_full(:,3:n_r)-h_H_12_alt(:,3:n_r))))))
        call writo('maximum error in h_H_33 = '//&
            &trim(r2str(maxval(abs(h_H_33_full(:,3:n_r)-h_H_33_alt(:,3:n_r))))))
        
        !! plot results
        !call print_GP_3D('h_H_11 (HEL,alt)','',reshape([h_H_11_full(:,3:n_r),&
            !&h_H_11_alt(:,3:n_r)],[nchi,n_r-2,2]),&
            !&x=reshape([R_H(:,3:n_r),R_H(:,3:n_r)],[nchi,n_r-2,2]),&
            !&y=reshape([Z_H(:,3:n_r),Z_H(:,3:n_r)],[nchi,n_r-2,2]))
        !call print_GP_3D('diff h_H_11','',h_H_11_full(:,3:n_r)-&
            !&h_H_11_alt(:,3:n_r),x=R_H(:,3:n_r),y=Z_H(:,3:n_r))
        !call print_GP_3D('h_H_12 (HEL,alt)','',reshape([h_H_12_full(:,3:n_r),&
            !&h_H_12_alt(:,3:n_r)],[nchi,n_r-2,2]),&
            !&x=reshape([R_H(:,3:n_r),R_H(:,3:n_r)],[nchi,n_r-2,2]),&
            !&y=reshape([Z_H(:,3:n_r),Z_H(:,3:n_r)],[nchi,n_r-2,2]))
        !call print_GP_3D('diff h_H_12','',h_H_12_full(:,3:n_r)-&
            !&h_H_12_alt(:,3:n_r),x=R_H(:,3:n_r),y=Z_H(:,3:n_r))
        !call print_GP_3D('h_H_33 (HEL,alt)','',reshape([h_H_33_full(:,3:n_r),&
            !&h_H_33_alt(:,3:n_r)],[nchi,n_r-2,2]),&
            !&x=reshape([R_H(:,3:n_r),R_H(:,3:n_r)],[nchi,n_r-2,2]),&
            !&y=reshape([Z_H(:,3:n_r),Z_H(:,3:n_r)],[nchi,n_r-2,2]))
        !call print_GP_3D('diff h_H_33','',h_H_33_full(:,3:n_r)-&
            !&h_H_33_alt(:,3:n_r),x=R_H(:,3:n_r),y=Z_H(:,3:n_r))
    end function check_metrics
    
    ! deallocates  HELENA  quantities  that  are not  used  any  more after  the
    ! equilibrium phase
    subroutine dealloc_HEL
        deallocate(h_H_11)
        deallocate(h_H_12)
        deallocate(h_H_33)
    end subroutine dealloc_HEL
    
    ! deallocates HELENA quantities that are not used any more
    subroutine dealloc_HEL_final
        deallocate(chi_H)
        deallocate(flux_H)
        deallocate(p0)
        deallocate(qs)
        deallocate(RBphi)
        deallocate(h_H_11_full)
        deallocate(h_H_12_full)
        deallocate(h_H_33_full)
        deallocate(R_H)
        deallocate(Z_H)
    end subroutine dealloc_HEL_final
end module HELENA
