!------------------------------------------------------------------------------!
!   Variables that have to do with HELENA quantities                           !
!------------------------------------------------------------------------------!
module HELENA_vars
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_1_type, X_2_type
    
    implicit none
    private
    public read_HEL, dealloc_HEL, &
        &pres_H, qs, flux_p_H, nchi, chi_H, ias, h_H_11, h_H_12, h_H_33, &
        &RBphi, R_H, Z_H
    
    ! global variables
    real(dp), allocatable :: chi_H(:)                                           ! poloidal angle
    real(dp), allocatable :: flux_p_H(:)                                        ! normal coordinate values (poloidal flux)
    real(dp), allocatable :: pres_H(:)                                          ! pressure profile
    real(dp), allocatable :: qs(:)                                              ! safety factor
    real(dp), allocatable :: RBphi(:)                                           ! R B_phi (= F)
    real(dp), allocatable :: h_H_11(:,:)                                        ! adapted upper metric factor 11 (gem11)
    real(dp), allocatable :: h_H_12(:,:)                                        ! adapted upper metric factor 12 (gem12)
    real(dp), allocatable :: h_H_33(:,:)                                        ! adapted upper metric factor 33 (1/gem33)
    real(dp), allocatable :: R_H(:,:)                                           ! major radius R (xout)
    real(dp), allocatable :: Z_H(:,:)                                           ! height Z (yout)
    integer :: nchi                                                             ! nr. of poloidal points (nchi)
    integer :: ias                                                              ! 0 if top-bottom symmetric, 1 if not

contains
    ! Reads the HELENA equilibrium data
    ! (from HELENA routine IODSK)
    ! Note: The variables in HELENA output are normalized globally w.r.t.
    !   - R_m: radius of magnetic axis at equilibrium
    !   - B_m: magnetic field at magnetic axis at equilibrium ,
    ! These have to be specified (see global variables above)
    ! Moreover, the variables R_H and Z_H  are also normalized w.r.t. a radius a
    ! and a diameter R_0:
    !   - R_H = (R-R_0)/a
    !   - Z_H = Z/a ,
    ! where R_0 and a are found through the variable "radius" and "eps"
    !   - radius = a / R_m
    !   - eps = a / R_0 ,
    ! as a function of R_m, which is the global normalization factor used here.
    ! Furthermore,  the  varaible  "cs"  contains the  sqrt  of  the  normalized
    ! flux/2pi  on  the normal  positions,  where  the normalization  factor  is
    ! contained in "cpsurf", which is the poloidal flux at the surface.
    ! Finally, some variables are not tabulated on the magnetic axis.
    ! Note that in the asymmetric case,  the number of poloidal points, nchi, is
    ! aumented by one  and the information for 0 is  copied into the information
    ! for 2pi, to facilitate interpolation.
    integer function read_HEL(n_r_eq,use_pol_flux_H) result(ierr)
        use num_vars, only: eq_name, eq_i
        
        character(*), parameter :: rout_name = 'read_HEL'
        
        ! input / output
        integer, intent(inout) :: n_r_eq                                        ! nr. of normal points in equilibrium grid
        logical, intent(inout) :: use_pol_flux_H                                ! .true. if HELENA equilibrium is based on pol. flux
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, kd                                                       ! counters
        integer :: nchi_loc                                                     ! local nchi
        real(dp), allocatable :: s_H(:)                                         ! flux coordinate s
        real(dp), allocatable :: dqs(:)                                         ! derivative of q
        real(dp), allocatable :: curj(:)                                        ! toroidal current
        real(dp), allocatable :: vx(:), vy(:)                                   ! R and Z of surface
        real(dp) :: Dj0, Dje                                                    ! derivative of toroidal current on axis and surface
        real(dp) :: Dpres_H_0, Dpres_H_e                                        ! normal derivative of pressure on axis and surface
        real(dp) :: dRBphi0, dRBphie                                            ! normal derivative of R B_phi on axis and surface
        real(dp) :: radius, raxis                                               ! minor radius, major radius
        real(dp) :: eps                                                         ! inverse aspect ratio
        real(dp) :: cpsurf                                                      ! poloidal flux on surface
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Reading data from HELENA output "' &
            &// trim(eq_name) // '"')
        call lvl_ud(1)
        
        ! set error messages
        err_msg = 'Failed to read HELENA output file'
        
        ! rewind input file
        rewind(eq_i)
        
        ! read mapped equilibrium from disk
        read(eq_i,*,IOSTAT=ierr) n_r_eq,ias                                     ! nr. normal points (JS0), top-bottom symmetry
        CHCKERR(err_msg)
        n_r_eq = n_r_eq + 1                                                     ! HELENA outputs nr. of normal points - 1
        
        allocate(s_H(n_r_eq))                                                   ! flux coordinate
        read(eq_i,*,IOSTAT=ierr) (s_H(kd),kd=1,n_r_eq)                          ! it is squared below, after reading cpsurf
        CHCKERR(err_msg)
        
        allocate(qs(n_r_eq))                                                    ! safety factor
        read(eq_i,*,IOSTAT=ierr) (qs(kd),kd=1,n_r_eq)
        CHCKERR(err_msg)
        
        allocate(dqs(n_r_eq))                                                   ! derivative of safety factor
        read(eq_i,*,IOSTAT=ierr) dqs(1),dqs(n_r_eq)                             ! first point, last point
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) (dqs(kd),kd=2,n_r_eq)                          ! second to last point (again)
        CHCKERR(err_msg)
        
        allocate(curj(n_r_eq))                                                  ! toroidal current
        read(eq_i,*,IOSTAT=ierr) (curj(kd),kd=1,n_r_eq)
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) Dj0,Dje                                        ! derivative of toroidal current at first and last point
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) nchi                                           ! nr. poloidal points
        CHCKERR(err_msg)
        nchi_loc = nchi
        if (ias.ne.0) nchi = nchi + 1                                           ! extend grid to 2pi
        
        allocate(chi_H(nchi))                                                   ! poloidal points
        read(eq_i,*,IOSTAT=ierr) (chi_H(id),id=1,nchi_loc)
        CHCKERR(err_msg)
        if (ias.ne.0) chi_H(nchi) = 2*pi
        
        allocate(h_H_11(nchi,n_r_eq))                                           ! upper metric factor 1,1
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_11(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (gem11)
        CHCKERR(err_msg)
        h_H_11(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        if (ias.ne.0) h_H_11(nchi,:) = h_H_11(1,:)
        
        allocate(h_H_12(nchi,n_r_eq))                                           ! upper metric factor 1,2
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_12(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (gem12)
        CHCKERR(err_msg)
        h_H_12(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        if (ias.ne.0) h_H_12(nchi,:) = h_H_12(1,:)
        
        read(eq_i,*,IOSTAT=ierr) cpsurf, radius                                 ! poloidal flux on surface, minor radius
        CHCKERR(err_msg)
        
        allocate(h_H_33(nchi,n_r_eq))                                           ! upper metric factor 3,3
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_33(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),id=&
            &nchi_loc+1,n_r_eq*nchi_loc)                                        ! (gem33)
        h_H_33(:,:) = 1._dp/h_H_33(:,:)                                         ! HELENA gives R^2, but need 1/R^2
        h_H_33(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        if (ias.ne.0) h_H_33(nchi,:) = h_H_33(1,:)
        
        read(eq_i,*,IOSTAT=ierr) raxis                                          ! major radius
        CHCKERR(err_msg)
        
        allocate(pres_H(n_r_eq))                                                ! pressure profile
        read(eq_i,*,IOSTAT=ierr) (pres_H(kd),kd=1,n_r_eq)
        CHCKERR(err_msg)
        read(eq_i,*,IOSTAT=ierr) Dpres_H_0,Dpres_H_e                            ! derivarives of pressure on axis and surface
        CHCKERR(err_msg)
        
        allocate(RBphi(n_r_eq))                                                 ! R B_phi (= F)
        read(eq_i,*,IOSTAT=ierr) (RBphi(kd),kd=1,n_r_eq)
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) dRBphi0,dRBphie                                ! derivatives of R B_phi on axis and surface
        CHCKERR(err_msg)
        
        allocate(vx(nchi))                                                      ! R B_phi
        read(eq_i,*,IOSTAT=ierr) (vx(id),id=1,nchi_loc)                         ! R on surface
        CHCKERR(err_msg)
        if (ias.ne.0) vx(nchi) = vx(1)
        
        allocate(vy(nchi))                                                      ! R B_phi
        read(eq_i,*,IOSTAT=ierr) (vy(id),id=1,nchi_loc)                         ! Z on surface
        CHCKERR(err_msg)
        if (ias.ne.0) vy(nchi) = vy(1)
        
        read(eq_i,*,IOSTAT=ierr) eps                                            ! inverse aspect ratio
        CHCKERR(err_msg)
        
        allocate(R_H(nchi,n_r_eq))                                              ! major radius R
        read(eq_i,*,IOSTAT=ierr) (R_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (xout)
        CHCKERR(err_msg)
        R_H(:,1) = R_H(:,2)                                                     ! first point is not given: set equal to second one
        if (ias.ne.0) R_H(nchi,:) = R_H(1,:)
        
        allocate(Z_H(nchi,n_r_eq))                                              ! height Z
        read(eq_i,*,IOSTAT=ierr) (Z_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (yout)
        CHCKERR(err_msg)
        Z_H(:,1) = Z_H(:,2)                                                     ! first point is not given: set equal to second one
        if (ias.ne.0) Z_H(nchi,:) = Z_H(1,:)
        
        ! transform back to unnormalized quantities
        radius = radius * raxis                                                 ! global length normalization with R_m
        flux_p_H = 2*pi*s_H**2 * cpsurf                                         ! rescale flux coordinate (HELENA uses psi_pol/2pi as flux)
        R_H = radius*(1._dp/eps + R_H)                                          ! local normalization with a
        Z_H = radius*Z_H                                                        ! local normalization with a
       
        ! HELENA always uses the poloidal flux
        use_pol_flux_H = .true.
        
        ! user output
        call writo('HELENA output given on '//trim(i2str(nchi))//&
            &' poloidal and '//trim(i2str(n_r_eq))//' normal points')
        call lvl_ud(-1)
        call writo('Data from HELENA output succesfully read')
        
        ! close the HELENA file
        close(eq_i)
    end function read_HEL
    
    ! deallocates HELENA quantities that are not used any more
    subroutine dealloc_HEL
        if (allocated(h_H_11)) deallocate(h_H_11)
        if (allocated(h_H_12)) deallocate(h_H_12)
        if (allocated(h_H_33)) deallocate(h_H_33)
        if (allocated(chi_H)) deallocate(chi_H)
        if (allocated(flux_p_H)) deallocate(flux_p_H)
        if (allocated(pres_H)) deallocate(pres_H)
        if (allocated(qs)) deallocate(qs)
        if (allocated(RBphi)) deallocate(RBphi)
        if (allocated(R_H)) deallocate(R_H)
        if (allocated(Z_H)) deallocate(Z_H)
    end subroutine dealloc_HEL
end module HELENA_vars
