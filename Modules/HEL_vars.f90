!------------------------------------------------------------------------------!
!   functions and routines that concern variables given by HELENA              !
!------------------------------------------------------------------------------!
module HEL_vars
#include <PB3D_macros.h>
    use message_ops, only: writo, lvl_ud
    use num_vars, only: dp, max_str_ln
    use output_ops, only: print_GP_2D, print_GP_3D
    
    implicit none
    private
    public read_HEL, dealloc_HEL, &
        &R_H, p0, qs, flux_H, nchi, chi, ias, h_H_11, h_H_12, h_H_33, RBphi
    
    ! global variables
    real(dp) :: R_H = 1.0_dp                                                    ! R of magnetic axis (normalization constant)
    real(dp) :: B_H = 1.0_dp                                                    ! B at magnetic axis (normalization constant)
    real(dp), allocatable :: p0(:)                                              ! pressure profile
    real(dp), allocatable :: qs(:)                                              ! safety factor
    real(dp), allocatable :: flux_H(:)                                          ! normal coordinate values
    real(dp), allocatable :: chi(:)                                             ! poloidal angle
    real(dp), allocatable :: RBphi(:)                                           ! R B_phi
    real(dp), allocatable :: h_H_11(:,:)                                        ! upper metric factor 11 (gem11)
    real(dp), allocatable :: h_H_12(:,:)                                        ! upper metric factor 12 (gem12)
    real(dp), allocatable :: h_H_33(:,:)                                        ! upper metric factor 33 (1/gem33)
    integer :: nchi                                                             ! nr. of poloidal points (nchi)
    integer :: ias                                                              ! 0 if top-bottom symmetric, 1 if not

contains
    ! Reads the HELENA equilibrium data
    ! (from HELENA routine IODSK)
    ! Note: The variables in HELENA output are normalized wrt.:
    !   - R_m: radius of magnetic axis at equilibrium
    !   - B_m: magnetic field at magnetic axis at equilibrium ,
    ! These have to be specified (see global variables above)
    ! However, the variables xout and yout are also normalized wrt.:
    !   - xout = (R-R_0)/a
    !   - R_0 yout Z/a ,
    ! where R_0 and a are found through the variable "radius" and "eps":
    !   - radius = a / R_m
    !   - eps = a / R_0 ,
    ! which link R_0 and a to R_m.
    ! Furthermore, the varaible "cs" contains the sqrt of the normalized flux on
    ! the  normal positions,  where  the normalization  factor  is contained  in
    ! "cpsurf", which is the poloidal flux at the surface.
    ! Finally, some variables are not tabulated on the magnetic axis.
    ! [MPI] only global master
    integer function read_HEL(n_r_eq,eq_use_pol_flux) result(ierr)
        use num_vars, only: glb_rank, eq_name, eq_i
        
        character(*), parameter :: rout_name = 'read_HEL'
        
        ! input / output
        integer, intent(inout) :: n_r_eq                                        ! nr. of normal points in equilibrium grid
        logical, intent(inout) :: eq_use_pol_flux                               ! .true. if equilibrium is based on poloidal flux
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        character(len=max_str_ln) :: err_msg_ias                                ! error message concerning ias
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: dqs(:)                                         ! derivative of q
        real(dp), allocatable :: curj(:)                                        ! toroidal current
        real(dp), allocatable :: xout(:,:)                                      ! major radius R
        real(dp), allocatable :: yout(:,:)                                      ! height Z
        real(dp), allocatable :: vx(:), vy(:)                                   ! R and Z of surface
        real(dp), allocatable :: jac_H(:,:)                                     ! Jacobian in (nchi,n_r_eq) variables
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
            err_msg_ias = 'You need to have a version of HELENA that outputs &
                &the variable IAS'
            
            ! rewind input file
            rewind(eq_i)
            
            ! read mapped equilibrium from disk
            read(eq_i,*,IOSTAT=ierr) n_r_eq,ias                                 ! nr. normal points (JS0), top-bottom symmetry
            CHCKERR(err_msg_ias)
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
            
            allocate(chi(nchi))                                                 ! poloidal points
            read(eq_i,*,IOSTAT=ierr) (chi(id),id=1,nchi)
            CHCKERR(err_msg)
            
            allocate(h_H_11(nchi,n_r_eq))
            read(eq_i,*,IOSTAT=ierr) (h_H_11(mod(id-1,nchi)+1,(id-1)/nchi+1),&
                &id=nchi+1,n_r_eq*nchi)                                         ! (gem11)
            CHCKERR(err_msg)
            h_H_11(:,1) = 0._dp                                                 ! first normal point is not given, so set to zero
            
            allocate(h_H_12(nchi,n_r_eq))
            read(eq_i,*,IOSTAT=ierr) (h_H_12(mod(id-1,nchi)+1,(id-1)/nchi+1),&
                &id=nchi+1,n_r_eq*nchi)                                         ! (gem12)
            CHCKERR(err_msg)
            h_H_12(:,1) = 0._dp                                                 ! first normal point is not given, so set to zero
            write(*,*) '!!!NOT SURE IF LHS -> RHS IS NECESSARY!!!'
            write(*,*) '!!!! IS THIS WITH FLUX OR WITH S ????????!??!?!??!'
            
            read(eq_i,*,IOSTAT=ierr) cpsurf, radius                             ! poloidal flux on surface, minor radius
            CHCKERR(err_msg)
            cpsurf = cpsurf * R_H**2 * B_H                                      ! back to real space
            flux_H = flux_H**2 * cpsurf                                         ! rescale poloidal flux
            
            allocate(h_H_33(nchi,n_r_eq))
            read(eq_i,*,IOSTAT=ierr) (h_H_33(mod(id-1,nchi)+1,(id-1)/nchi+1),&
                &id=nchi+1,n_r_eq*nchi)                                         ! (gem11)
            h_H_33(:,:) = 1/h_H_33(:,:)                                         ! HELENA gives R^2, but need 1/R^2
            h_H_33(:,1) = 0._dp                                                 ! first normal point is not given, so set to zero
            
            read(eq_i,*,IOSTAT=ierr) raxis                                      ! major radius
            CHCKERR(err_msg)
            
            allocate(p0(n_r_eq))                                                ! pressure profile
            read(eq_i,*,IOSTAT=ierr) (p0(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            
            read(eq_i,*,IOSTAT=ierr) dp0,dpe                                    ! derivarives of pressure on axis and surface
            CHCKERR(err_msg)
            
            allocate(RBphi(n_r_eq))                                             ! R B_phi
            read(eq_i,*,IOSTAT=ierr) (RBphi(kd),kd=1,n_r_eq)
            CHCKERR(err_msg)
            
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
            
            allocate(xout(nchi,n_r_eq))                                         ! major radius R
            xout(:,1) = 0._dp                                                   ! values on axis are not given by HELENA -> set to zero
            read(eq_i,*,IOSTAT=ierr) &
                &(xout(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)   ! (gem11)
            CHCKERR(err_msg)
            xout = R_H*radius*(1._dp/eps + xout)                                ! back to real space
            
            allocate(yout(nchi,n_r_eq))                                         ! height Z
            yout(:,1) = 0._dp                                                   ! values on axis are not given by HELENA -> set to zero
            read(eq_i,*,IOSTAT=ierr) &
                &(yout(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)   ! (gem11)
            CHCKERR(err_msg)
            yout = R_H*radius*yout                                              ! back to real space
            
            ! set up Jacobian
            allocate(jac_H(nchi,n_r_eq))
            do kd = 1,n_r_eq
                jac_H(:,kd) = 1.0_dp/(h_H_33(:,kd)*qs(kd)*RBphi(kd))
            end do
            
            write(*,*) 'GO BACK TO REAL SPACE !!!!!!!!!!!!!1'
            
            ! HELENA always uses the poloidal flux
            eq_use_pol_flux = .true.
            
            call lvl_ud(-1)
            call writo('Grid parameters successfully read')
        end if
    end function read_HEL
    
    ! deallocates HELENA quantities that are not used anymore
    subroutine dealloc_HEL
    end subroutine dealloc_HEL
end module HEL_vars
