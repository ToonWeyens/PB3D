!------------------------------------------------------------------------------!
!   operations and variable pertaining to the HELENA equilibrim                !
!------------------------------------------------------------------------------!
module HELENA
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use met_vars, only: met_type
    use X_vars, only: X_type
    
    implicit none
    private
    public read_HEL, dealloc_HEL, interp_HEL_on_grid, &
        &pres_H, qs, flux_p_H, nchi, chi_H, ias, h_H_11, h_H_12, h_H_33, &
        &RBphi, R_H, Z_H, R_0_H, B_0_H
#if ldebug
    public test_metrics_H
#endif
    
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
    real(dp), parameter :: R_0_H = 2.0_dp                                       ! normalization radius
    real(dp), parameter :: B_0_H = 2.0_dp                                       ! normalization magnetic field
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
    ! [MPI] only global master
    integer function read_HEL(n_r_eq,use_pol_flux_H) result(ierr)
        use num_vars, only: eq_name, eq_i
        
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
        
        allocate(flux_p_H(n_r_eq))                                              ! flux, derived from normal coordinate
        read(eq_i,*,IOSTAT=ierr) (flux_p_H(kd),kd=1,n_r_eq)                     ! it is squared below, after reading cpsurf
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
        
        allocate(chi_H(nchi))                                                   ! poloidal points
        read(eq_i,*,IOSTAT=ierr) (chi_H(id),id=1,nchi)
        CHCKERR(err_msg)
        
        allocate(h_H_11(nchi,n_r_eq))                                           ! upper metric factor 1,1
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_11(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,&
            &n_r_eq*nchi)                                                       ! (gem11)
        CHCKERR(err_msg)
        h_H_11(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        
        allocate(h_H_12(nchi,n_r_eq))                                           ! upper metric factor 1,2
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_12(mod(id-1,nchi)+1,(id-1)/nchi+1),&
            &id=nchi+1,n_r_eq*nchi)                                             ! (gem12)
        CHCKERR(err_msg)
        h_H_12(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        
        read(eq_i,*,IOSTAT=ierr) cpsurf, radius                                 ! poloidal flux on surface, minor radius
        CHCKERR(err_msg)
        flux_p_H = 2*pi*flux_p_H**2 * cpsurf                                    ! rescale poloidal flux (HELENA uses psi_pol/2pi as flux)
        
        allocate(h_H_33(nchi,n_r_eq))                                           ! upper metric factor 3,3
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_33(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,&
            &n_r_eq*nchi)                                                       ! (gem33)
        h_H_33(:,:) = 1._dp/h_H_33(:,:)                                         ! HELENA gives R^2, but need 1/R^2
        h_H_33(:,1) = 0._dp                                                     ! first normal point is not given, so set to zero
        
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
        read(eq_i,*,IOSTAT=ierr) (vx(id),id=1,nchi)                             ! R on surface
        CHCKERR(err_msg)
        
        allocate(vy(nchi))                                                      ! R B_phi
        read(eq_i,*,IOSTAT=ierr) (vy(id),id=1,nchi)                             ! Z on surface
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) eps                                            ! inerse aspect ratio
        CHCKERR(err_msg)
        
        allocate(R_H(nchi,n_r_eq))                                              ! major radius R
        R_H(:,1) = 0._dp                                                        ! values on axis are not given by HELENA -> set to zero
        read(eq_i,*,IOSTAT=ierr) &
            &(R_H(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)        ! (xout)
        CHCKERR(err_msg)
        R_H = radius*(1._dp/eps + R_H)                                          ! from a to R
        
        allocate(Z_H(nchi,n_r_eq))                                              ! height Z
        Z_H(:,1) = 0._dp                                                        ! values on axis are not given by HELENA -> set to zero
        read(eq_i,*,IOSTAT=ierr) &
            &(Z_H(mod(id-1,nchi)+1,(id-1)/nchi+1),id=nchi+1,n_r_eq*nchi)        ! (yout)
        CHCKERR(err_msg)
        Z_H = radius*Z_H                                                        ! from a to R
       
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
        deallocate(h_H_11)
        deallocate(h_H_12)
        deallocate(h_H_33)
        deallocate(chi_H)
        deallocate(flux_p_H)
        deallocate(pres_H)
        deallocate(qs)
        deallocate(RBphi)
        deallocate(R_H)
        deallocate(Z_H)
    end subroutine dealloc_HEL
    
    ! calculate interpolation  factors for angular interpolation  in grid_out of
    ! quantities defined on grid_in. This version  is specific for an input grid
    ! corresponding to  axisymmetric variables with optional  top-down symmetry,
    ! as is the case for variables resulting from HELENA equilibria.
    ! The output of a 3D array of real values for the poloidal angle theta where
    ! the  floored  integer of  each  value  indicates  the  base index  of  the
    ! interpolated value in the output grid  and the modulus is the the fraction
    ! towards the next integer.
    ! The flag  td_sym indicates optionally  that there is top-down  symmetry as
    ! well as axisymmetry. When there is top-down symmetry, the variables in the
    ! lower half  (i.e. pi<theta<2pi) are  calculated from the variables  in the
    ! upper  half  using  var(2pi-theta)  = var(theta).  However,  some  complex
    ! variables also  need to  apply the complex  conjugate when  using top-down
    ! symmetry.  Therefore, to  indicate  this, the  sign  of the  interpolation
    ! factor is inverted so  it is negative. For variables that  do not need the
    ! complex conjugate, this can safely be ignored by taking the absolute value
    ! of the interpolation factor.
    ! By default, the variables in the Flux coord. system are used, but this can
    ! be changed optionally with the flag "use_E".
    integer function get_ang_interp_data_HEL(grid_in,grid_out,theta_i,td_sym,&
        &use_E) result(ierr)
        use utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'get_ang_interp_data_HEL'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in, grid_out                        ! input and output grid
        real(dp), allocatable, intent(inout) :: theta_i(:,:,:)                  ! interpolation index
        logical, intent(in), optional :: td_sym                                 ! top-down symmetry
        logical, intent(in), optional :: use_E                                  ! whether E is used instead of F
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        logical :: td_sym_loc                                                   ! local version of td_sym
        logical :: use_E_loc                                                    ! local version of use_E
        real(dp) :: theta_loc                                                   ! local theta of output grid
        real(dp), pointer :: theta_in(:) => null()                              ! input theta
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! test whether axisymmetric grid
        if (grid_in%n(2).ne.1) then
            ierr = 1
            err_msg = 'Not an axisymmetric grid'
            CHCKERR(err_msg)
        end if
        ! test whether normal sizes compatible
        if (grid_in%grp_n_r.ne.grid_out%grp_n_r) then
            ierr = 1
            err_msg = 'Grids are not compatible in normal direction'
            CHCKERR(err_msg)
        end if
        
        ! set up local use_E and td_sym
        use_E_loc = .false.
        if (present(use_E)) use_E_loc = use_E
        td_sym_loc = .false.
        if (present(td_sym)) td_sym_loc = td_sym
        
        ! allocate theta_i
        allocate(theta_i(grid_out%n(1),grid_out%n(2),grid_out%grp_n_r))
        
        ! For every point on output grid, check to which half poloidal circle on
        ! the  axisymmetric input  grid  it  belongs to.  If  there is  top-down
        ! symmetry and the  angle theta lies in the bottom  part, the quantities
        ! have to be taken from their symmetric counterpart (2pi-theta).
        ! loop over all normal points of this rank
        do kd = 1,grid_out%grp_n_r
            ! point theta_in
            if (use_E) then
                theta_in => grid_in%theta_E(:,1,kd)                             ! axisymmetric grid should not have only one toroidal point
            else
                theta_in => grid_in%theta_F(:,1,kd)                             ! axisymmetric grid should not have only one toroidal point
            end if
            
            ! loop over all angles of ang_2
            do jd = 1,grid_out%n(2)
                ! loop over all angles of ang_1
                do id = 1,grid_out%n(1)
                    ! set the local poloidal point from theta
                    if (use_E) then
                        theta_loc = grid_out%theta_E(id,jd,kd)
                    else
                        theta_loc = grid_out%theta_F(id,jd,kd)
                    end if
                    
                    ! add or subtract  2pi to the parallel angle until  it is at
                    ! least 0 to get principal range 0..2pi
                    if (theta_loc.lt.0._dp) then
                        do while (theta_loc.lt.0._dp)
                            theta_loc = theta_loc + 2*pi
                        end do
                    else if (theta_loc.gt.2*pi) then
                        do while (theta_loc.gt.2._dp*pi)
                            theta_loc = theta_loc - 2*pi
                        end do
                    end if
                    
                    ! get interpolation factors in theta_i
                    if (td_sym .and. theta_loc.gt.pi) then                      ! bottom part and top-down symmetric
                        ierr = con2dis(2*pi-theta_loc,theta_i(id,jd,kd),&
                            &theta_in)
                        theta_i(id,jd,kd) = - theta_i(id,jd,kd)                 ! top-down symmetry applied
                    else
                        ierr = con2dis(theta_loc,theta_i(id,jd,kd),theta_in)
                    end if
                    CHCKERR('')
                end do
            end do
        end do
        
        ! clean up
        nullify(theta_in)
    end function get_ang_interp_data_HEL
    
    ! Adapt  some variables  resulting from  HELENA equilibria  to field-aligned
    ! grid (angularly):
    !   - perturbation variables: U_i, DU_i, PV_i, KV_i
    ! and optionally:
    !   - equilibrium variables of interest are flux variables
    !   - metric variables: jac_FD, g_FD, h_FD
    integer function interp_HEL_on_grid(grid_eq,grid_eq_B,X,X_B,met,met_B,&
        &eq,eq_B,grid_name) result(ierr)
        use num_vars, only: prog_style
        
        character(*), parameter :: rout_name = 'interp_HEL_on_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq, grid_eq_B                       ! general and field-aligned equilibrium grid
        type(X_type), intent(in) :: X                                           ! general perturbation variables
        type(X_type), intent(inout) :: X_B                                      ! field-aligned perturbation variables
        type(met_type), intent(in), optional :: met                             ! general metric variables
        type(met_type), intent(inout), optional :: met_B                        ! field-aligned metric variables
        type(eq_type), intent(in), optional :: eq                               ! general equilibrium variables
        type(eq_type), intent(inout), optional :: eq_B                          ! field-aligned equilibrium variables
        character(len=*), intent(in), optional :: grid_name                     ! name of grid to which to adapt quantities
        
        ! local variables
        real(dp), allocatable :: theta_i(:,:,:)                                 ! interpolation index
        logical :: td_sym                                                       ! whether there is top-down symmetry (relevant only for HELENA)
        character(len=max_str_ln) :: grid_name_loc                              ! local version of grid_name
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (prog_style.eq.2) then                                               ! PB3D_POST
            if (.not.present(eq) .or. .not.present(eq_B)) then
                ierr = 1
                err_msg = 'For PB3D_POST, eq and eq_B are needed as well'
                CHCKERR(err_msg)
            end if
        end if
        
        ! set up local grid name
        grid_name_loc = ''
        if (present(grid_name)) grid_name_loc = ' to '//grid_name
        
        ! user output
        call writo('Adapting quantities'//trim(grid_name_loc))
        call lvl_ud(1)
        
        ! set up td_sym
        if (ias.eq.0) then
            td_sym = .true.
        else
            td_sym = .false.
        end if
        
        ! get angular interpolation factors
        ierr = get_ang_interp_data_HEL(grid_eq,grid_eq_B,theta_i,td_sym=td_sym,&
            &use_E=.false.)
        CHCKERR('')
        
        ! user output
        call writo('Adapting perturbation quantities')
        call lvl_ud(1)
        ! adapt common variables for all program styles
        X_B%vac_res = X%vac_res
        call interp_var_4D_complex(X%U_0,theta_i,X_B%U_0,sym_type=2)
        call interp_var_4D_complex(X%U_1,theta_i,X_B%U_1,sym_type=2)
        call interp_var_4D_complex(X%DU_0,theta_i,X_B%DU_0,sym_type=1)
        call interp_var_4D_complex(X%DU_1,theta_i,X_B%DU_1,sym_type=1)
        ! adapt custom variables depending on program style
        select case (prog_style)
            case(1)                                                             ! PB3D
                call interp_var_4D_complex(X%J_exp_ang_par_F,theta_i,&
                    &X_B%J_exp_ang_par_F,sym_type=1)
                call interp_var_4D_complex(X%PV_0,theta_i,X_B%PV_0,sym_type=1)
                call interp_var_4D_complex(X%PV_1,theta_i,X_B%PV_1,sym_type=1)
                call interp_var_4D_complex(X%PV_2,theta_i,X_B%PV_2,sym_type=1)
                call interp_var_4D_complex(X%KV_0,theta_i,X_B%KV_0,sym_type=1)
                call interp_var_4D_complex(X%KV_1,theta_i,X_B%KV_1,sym_type=1)
                call interp_var_4D_complex(X%KV_2,theta_i,X_B%KV_2,sym_type=1)
            case(2)                                                             ! PB3D_POST
                ! do nothing
            case default
                err_msg = 'No program style associated with '//&
                    &trim(i2str(prog_style))
                ierr = 1
                CHCKERR(err_msg)
        end select
        call lvl_ud(-1)
        
        if (present(eq).and.present(eq_B)) then
            ! user output
            call writo('Adapting equilibrium quantities')
            call lvl_ud(1)
            eq_B%pres_FD = eq%pres_FD
            eq_B%q_saf_FD = eq%q_saf_FD
            eq_B%rot_t_FD = eq%rot_t_FD
            eq_B%flux_p_FD = eq%flux_p_FD
            eq_B%flux_t_FD = eq%flux_t_FD
            eq_B%rho = eq%rho
            call interp_var_3D_real(eq%S,theta_i,eq_B%S)
            call interp_var_3D_real(eq%sigma,theta_i,eq_B%sigma)
            call interp_var_3D_real(eq%kappa_n,theta_i,eq_B%kappa_n)
            call interp_var_3D_real(eq%kappa_g,theta_i,eq_B%kappa_g,&
                &sym_type=2)
            call lvl_ud(-1)
        end if
        
        if (present(met).and.present(met_B)) then
            ! user output
            call writo('Adapting metric quantities')
            call lvl_ud(1)
            call interp_var_6D_real(met%jac_FD,theta_i,met_B%jac_FD)
            call interp_var_7D_real(met%g_FD,theta_i,met_B%g_FD)
            call interp_var_7D_real(met%h_FD,theta_i,met_B%h_FD)
            call lvl_ud(-1)
        end if
        
        ! user output
        call lvl_ud(-1)
        call writo('Quantities adapted'//trim(grid_name_loc))
    contains
        ! Interpolate a  variable defined  on an  axisymmetric grid  at poloidal
        ! angles indicated by the interpolation factors theta_i.
        ! There  is  an  optional  variable   sym_type  that  allows  for  extra
        ! operations to be done on the variable if top-down symmetry is applied:
        !   - sym_type = 0: var(2pi-theta) = var(theta)
        !   - sym_type = 1: var(2pi-theta) = var(theta)*    (only for complex)
        !   - sym_type = 2: var(2pi-theta) = -var(theta)*
        ! When top-down  symmetry has been  used to calculate  the interpolation
        ! factors, this is indicated by a  negative factor instead of a positive
        ! one.
        ! (also  correct  if i_hi = i_lo)
        subroutine interp_var_3D_real(var,theta_i,var_int,sym_type)             ! 3D_real version
            ! input / output
            real(dp), intent(in) :: var(:,:,:)                                  ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:)                           ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd                                               ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! iterate over all normal points
            do kd = 1,size(var_int,3)
                ! iterate over all geodesical points
                do jd = 1,size(var_int,2)
                    ! iterate over all parallel points
                    do id = 1,size(var_int,1)
                        ! set up i_lo and i_hi
                        i_lo = floor(abs(theta_i(id,jd,kd)))
                        i_hi = ceiling(abs(theta_i(id,jd,kd)))
                        
                        var_int(id,jd,kd) = var(i_lo,1,kd) + &
                            &(var(i_hi,1,kd)-var(i_lo,1,kd))*&
                            &(abs(theta_i(id,jd,kd))-i_lo)                      ! because i_hi - i_lo = 1
                        if (theta_i(id,jd,kd).lt.0) then
                            if (sym_type_loc.eq.2) var_int(id,jd,kd) = &
                                &- var_int(id,jd,kd)
                        end if
                    end do
                end do
            end do
        end subroutine interp_var_3D_real
        subroutine interp_var_6D_real(var,theta_i,var_int,sym_type)             ! 6D_real version
            ! input / output
            real(dp), intent(in) :: var(:,:,:,:,:,:)                            ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:,:,:,:)                     ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd, ld                                           ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! iterate over all normal points
            do kd = 1,size(var_int,3)
                ! iterate over all geodesical points
                do jd = 1,size(var_int,2)
                    ! iterate over all parallel points
                    do id = 1,size(var_int,1)
                        ! set up i_lo and i_hi
                        i_lo = floor(abs(theta_i(id,jd,kd)))
                        i_hi = ceiling(abs(theta_i(id,jd,kd)))
                        
                        var_int(id,jd,kd,:,:,:) = var(i_lo,1,kd,:,:,:) + &
                            &(var(i_hi,1,kd,:,:,:)-var(i_lo,1,kd,:,:,:))*&
                            &(abs(theta_i(id,jd,kd))-i_lo)                      ! because i_hi - i_lo = 1
                        if (theta_i(id,jd,kd).lt.0) then
                            if (sym_type_loc.eq.2) var_int(id,jd,kd,:,:,:) = &
                                &- var_int(id,jd,kd,:,:,:)
                        end if
                    end do
                end do
            end do
            
            ! transform d/dtheta
            do ld = 1,size(var_int,6)
                var_int(:,:,:,:,:,ld) = (-1)**(ld-1)*var_int(:,:,:,:,:,ld)
            end do
        end subroutine interp_var_6D_real
        subroutine interp_var_7D_real(var,theta_i,var_int,sym_type)             ! 7D_real version
            ! input / output
            real(dp), intent(in) :: var(:,:,:,:,:,:,:)                          ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:,:,:,:,:)                   ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd, ld                                           ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! iterate over all normal points
            do kd = 1,size(var_int,3)
                ! iterate over all geodesical points
                do jd = 1,size(var_int,2)
                    ! iterate over all parallel points
                    do id = 1,size(var_int,1)
                        ! set up i_lo and i_hi
                        i_lo = floor(abs(theta_i(id,jd,kd)))
                        i_hi = ceiling(abs(theta_i(id,jd,kd)))
                        
                        var_int(id,jd,kd,:,:,:,:) = var(i_lo,1,kd,:,:,:,:) + &
                            &(var(i_hi,1,kd,:,:,:,:)-var(i_lo,1,kd,:,:,:,:))*&
                            &(abs(theta_i(id,jd,kd))-i_lo)                      ! because i_hi - i_lo = 1
                        if (theta_i(id,jd,kd).lt.0) then
                            if (sym_type_loc.eq.2) var_int(id,jd,kd,:,:,:,:) = &
                                &- var_int(id,jd,kd,:,:,:,:)
                        end if
                    end do
                end do
            end do
            
            ! transform d/dtheta
            do ld = 1,size(var_int,7)
                var_int(:,:,:,:,:,:,ld) = (-1)**(ld-1)*var_int(:,:,:,:,:,:,ld)
            end do
        end subroutine interp_var_7D_real
        subroutine interp_var_4D_complex(var,theta_i,var_int,sym_type)          ! 4D_complex version
            ! input / output
            complex(dp), intent(in) :: var(:,:,:,:)                             ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            complex(dp), intent(inout) :: var_int(:,:,:,:)                      ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd                                               ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! iterate over all normal points
            do kd = 1,size(var_int,3)
                ! iterate over all geodesical points
                do jd = 1,size(var_int,2)
                    ! iterate over all parallel points
                    do id = 1,size(var_int,1)
                        ! set up i_lo and i_hi
                        i_lo = floor(abs(theta_i(id,jd,kd)))
                        i_hi = ceiling(abs(theta_i(id,jd,kd)))
                        
                        var_int(id,jd,kd,:) = var(i_lo,1,kd,:) + &
                            &(var(i_hi,1,kd,:)-var(i_lo,1,kd,:))*&
                            &(abs(theta_i(id,jd,kd))-i_lo)                      ! because i_hi - i_lo = 1
                        if (theta_i(id,jd,kd).lt.0) then
                            if (sym_type_loc.eq.1) var_int(id,jd,kd,:) = &
                                &conjg(var_int(id,jd,kd,:))
                            if (sym_type_loc.eq.2) var_int(id,jd,kd,:) = &
                                &- conjg(var_int(id,jd,kd,:))
                        end if
                    end do
                end do
            end do
        end subroutine interp_var_4D_complex
    end function interp_HEL_on_grid
    
#if ldebug
    ! Checks whether the metric elements  provided by HELENA are consistent with
    ! a direct calculation using the coordinate transformations:
    !   |nabla psi|^2           = 1/jac^2 ((dZ/dchi)^2 + (dR/dchi)^2)
    !   |nabla psi nabla chi|   = 1/jac^2 (dZ/dchi dZ/dpsi + dR/dchi dR/dpsi)
    !   |nabla chi|^2           = 1/jac^2 ((dZ/dpsi)^2 + (dR/dpsi)^2)
    !   |nabla phi|^2           = 1/R^2
    ! with jac = dZ/dpsi dR/dchi - dR/dpsi dZ/dchi
    ! Also, test whether the pressure balance is satisfied.
    integer function test_metrics_H(n_r) result(ierr)
        use num_vars, only: glb_rank
        use utilities, only: calc_deriv
        use output_ops, only: plot_diff_HDF5
        
        character(*), parameter :: rout_name = 'test_metrics_H'
        
        ! input / output
        integer, intent(inout) :: n_r                                           ! nr. of normal points in equilibrium grid
        
        ! local variables
        integer :: id, kd                                                       ! counters
        real(dp), allocatable :: Rchi(:,:), Rpsi(:,:)                           ! chi and psi derivatives of R
        real(dp), allocatable :: Zchi(:,:), Zpsi(:,:)                           ! chi and psi derivatives of Z
        real(dp), allocatable :: jac(:,:)                                       ! jac as defined above
        real(dp), allocatable :: h_H_11_alt(:,:,:)                              ! alternative calculation for upper metric factor 11
        real(dp), allocatable :: h_H_12_alt(:,:,:)                              ! alternative calculation for upper metric factor 12
        real(dp), allocatable :: h_H_33_alt(:,:,:)                              ! alternative calculation for upper metric factor 33
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: r_min = 4                                                    ! first normal index that has meaning
        real(dp), allocatable :: tempvar(:,:,:,:)                               ! temporary variable
        
        ! initialize ierr
        ierr = 0
        
        if (glb_rank.eq.0) then
            ! user output
            call writo('Checking consistency of metric factors')
            call lvl_ud(1)
            
            ! calculate  the  auxiliary quantities  Zchi,  zpsi,  Rchi and  Rpsi
            ! containing the derivatives as well as jac
            allocate(Rchi(nchi,n_r),Rpsi(nchi,n_r))
            allocate(Zchi(nchi,n_r),Zpsi(nchi,n_r))
            allocate(jac(nchi,n_r))
            
            do id = 1,nchi
                ierr = calc_deriv(R_H(id,:),Rpsi(id,:),flux_p_H/(2*pi),1,2)
                CHCKERR('')
                ierr = calc_deriv(Z_H(id,:),Zpsi(id,:),flux_p_H/(2*pi),1,2)
                CHCKERR('')
            end do
            do kd = 1,n_r
                ierr = calc_deriv(R_H(:,kd),Rchi(:,kd),chi_H,1,2)
                CHCKERR('')
                ierr = calc_deriv(Z_H(:,kd),Zchi(:,kd),chi_H,1,2)
                CHCKERR('')
            end do
            jac = Zpsi*Rchi - Zchi*Rpsi
            
            ! calculate the metric factors directly
            allocate(h_H_11_alt(nchi,1,n_r))
            allocate(h_H_12_alt(nchi,1,n_r))
            allocate(h_H_33_alt(nchi,1,n_r))
            
            h_H_11_alt(:,1,:) = 1._dp/(jac**2) * (Zchi**2 + Rchi**2)
            h_H_12_alt(:,1,:) = -1._dp/(jac**2) * (Zchi*Zpsi + Rchi*Rpsi)
            h_H_33_alt(:,1,:) = 1._dp/(R_H**2)
            
            ! output h_H_11
            ! set some variables
            file_name = 'TEST_h_H_1_1'
            description = 'Testing whether the HELENA metric factor h_1_1 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_11_alt(:,:,r_min:n_r),&
                &reshape(h_H_11(:,r_min:n_r),[nchi,1,n_r-r_min+1]),&
                &file_name,description=description,output_message=.true.)
            
            ! output h_H_12
            ! set some variables
            file_name = 'TEST_h_H_1_2'
            description = 'Testing whether the HELENA metric factor h_1_2 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_12_alt(:,:,r_min:n_r),&
                &reshape(h_H_12(:,r_min:n_r),[nchi,1,n_r-r_min+1]),&
                &file_name,description=description,output_message=.true.)
            
            ! output h_H_33
            ! set some variables
            file_name = 'TEST_h_H_3_3'
            description = 'Testing whether the HELENA metric factor h_3_3 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_33_alt(:,:,r_min:n_r),&
                &reshape(h_H_33(:,r_min:n_r),[nchi,1,n_r-r_min+1]),&
                &file_name,description=description,output_message=.true.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
            
            ! user output
            call writo('Checking pressure balance')
            call lvl_ud(1)
            
            ! calculate auxiliary  quantities:
            !   1: D1 F ,
            !   2: D1 p ,
            !   3: D1 (q/F h_11) ,
            !   4: D2 (q/F h_12) .
            allocate(tempvar(nchi,1,n_r,4))
            do id = 1,nchi
                ierr = calc_deriv(RBphi,tempvar(id,1,:,1),flux_p_H/(2*pi),1,1)
                CHCKERR('')
                ierr = calc_deriv(pres_H,tempvar(id,1,:,2),flux_p_H/(2*pi),1,1)
                CHCKERR('')
                ierr = calc_deriv(qs/RBphi*h_H_11(id,:),tempvar(id,1,:,3),&
                    &flux_p_H/(2*pi),1,1)
                CHCKERR('')
            end do
            do kd = 1,n_r
                ierr = calc_deriv(qs(kd)/RBphi(kd)*h_H_12(:,kd),&
                    &tempvar(:,1,kd,4),chi_H,1,1)
                CHCKERR('')
            end do
            
            ! calculate pressure  balance in tempvar(1)
            !   mu_0 p' = F/(qR^2) (d/d1 (h_11 q/F) + d/d2 (h_12 q/F) + q F')
            do kd = 1,n_r
                tempvar(:,1,kd,1) = -RBphi(kd)*h_H_33(:,kd)/qs(kd) * &
                    &(tempvar(:,1,kd,1)*qs(kd) + tempvar(:,1,kd,3) + &
                    &tempvar(:,1,kd,4))
            end do
            
            ! output difference with p'
            ! set some variables
            file_name = 'TEST_p_H'
            description = 'Testing whether the HELENA pressure balance is &
                &consistent'
            
            ! plot difference
            call plot_diff_HDF5(tempvar(:,:,:,1),tempvar(:,:,:,2),&
                &file_name,description=description,output_message=.true.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end if
    end function test_metrics_H
#endif
end module HELENA
