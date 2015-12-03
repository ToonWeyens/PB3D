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
    use X_vars, only: X_1_type, X_2_type
    
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
        flux_p_H = 2*pi*s_H**2 * cpsurf                                         ! rescale flux coordinate (HELENA uses psi_pol/2pi as flux)
        
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
        
        read(eq_i,*,IOSTAT=ierr) eps                                            ! inerse aspect ratio
        CHCKERR(err_msg)
        
        allocate(R_H(nchi,n_r_eq))                                              ! major radius R
        R_H(:,1) = 0._dp                                                        ! values on axis are not given by HELENA -> set to zero
        read(eq_i,*,IOSTAT=ierr) (R_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (xout)
        CHCKERR(err_msg)
        R_H = radius*(1._dp/eps + R_H)                                          ! from a to R
        if (ias.ne.0) R_H(nchi,:) = R_H(1,:)
        
        allocate(Z_H(nchi,n_r_eq))                                              ! height Z
        Z_H(:,1) = 0._dp                                                        ! values on axis are not given by HELENA -> set to zero
        read(eq_i,*,IOSTAT=ierr) (Z_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_eq*nchi_loc)                                     ! (yout)
        CHCKERR(err_msg)
        Z_H = radius*Z_H                                                        ! from a to R
        if (ias.ne.0) Z_H(nchi,:) = Z_H(1,:)
       
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
    
    ! calculate interpolation  factors for angular interpolation  in grid_out of
    ! quantities defined on grid_in. This version  is specific for an input grid
    ! corresponding to  axisymmetric variables with optional  top-down symmetry,
    ! as is the case for variables resulting from HELENA equilibria.
    ! The output of a 3D array of real values for the poloidal angle theta where
    ! the  floored  integer of  each  value  indicates  the  base index  of  the
    ! interpolated value in the output grid  and the modulus is the the fraction
    ! towards the next integer.
    ! The flag tb_sym indicates optionally  that there is top-bottom symmetry as
    ! well as axisymmetry. When there is top-down symmetry, the variables in the
    ! lower half  (i.e. -pi<theta<0)  are calculated from  the variables  in the
    ! upper half  using the  symmetry properties of  the variables.  To indicate
    ! this,  the sign  of the  interpolation factor  is inverted  to a  negative
    ! value.
    ! The displacement of the fundamental theta interval is also outputted. For
    ! asymmetric variables this is for example:
    !    1  for -2pi ..  0
    !    0  for  0   ..  2pi
    !   -1  for  2pi ..  4pi
    ! etc.
    ! For symmetric variables, this is for example:
    !    1  for -3pi .. -2pi,   with symmetry property
    !    1  for -2pi .. -pi
    !    0  for -pi  ..  0,     with symmetry property
    !    0  for  0   ..  pi
    !   -1  for  pi  ..  2pi,   with symmetry property
    !   -1  for  3pi ..  3pi
    ! etc.
    ! By default, the variables in the Flux coord. system are used, but this can
    ! be changed optionally with the flag "use_E".
    integer function get_ang_interp_data_HEL(grid_in,grid_out,theta_i,&
        &fund_theta_int_displ,tb_sym,use_E) result(ierr)
        use utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'get_ang_interp_data_HEL'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in, grid_out                        ! input and output grid
        real(dp), allocatable, intent(inout) :: theta_i(:,:,:)                  ! interpolation index
        integer, allocatable, intent(inout) :: fund_theta_int_displ(:,:,:)      ! displacement of fundamental theta interval
        logical, intent(in), optional :: tb_sym                                 ! top-bottom symmetry
        logical, intent(in), optional :: use_E                                  ! whether E is used instead of F
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        logical :: tb_sym_loc                                                   ! local version of tb_sym
        logical :: use_E_loc                                                    ! local version of use_E
        real(dp) :: theta_loc                                                   ! local theta of output grid
        real(dp), pointer :: theta_in(:) => null()                              ! input theta
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: fund_theta_int_lo, fund_theta_int_hi                         ! lower and upper limit of fundamental theta interval
        
        ! initialize ierr
        ierr = 0
        
        ! test whether axisymmetric grid
        if (grid_in%n(2).ne.1) then
            ierr = 1
            err_msg = 'Not an axisymmetric grid'
            CHCKERR(err_msg)
        end if
        ! test whether normal sizes compatible
        if (grid_in%loc_n_r.ne.grid_out%loc_n_r) then
            ierr = 1
            err_msg = 'Grids are not compatible in normal direction'
            CHCKERR(err_msg)
        end if
        
        ! set up local use_E and tb_sym
        use_E_loc = .false.
        if (present(use_E)) use_E_loc = use_E
        tb_sym_loc = .false.
        if (present(tb_sym)) tb_sym_loc = tb_sym
        
        ! allocate theta_i and fundamental theta interval displacement
        allocate(theta_i(grid_out%n(1),grid_out%n(2),grid_out%loc_n_r))
        allocate(fund_theta_int_displ(grid_out%n(1),grid_out%n(2),&
            &grid_out%loc_n_r))
        
        ! initialize fundamental theta interval displacement
        fund_theta_int_displ = 0
        
        ! For every point on output grid, check to which half poloidal circle on
        ! the  axisymmetric input  grid  it  belongs to.  If  there is  top-down
        ! symmetry and the  angle theta lies in the bottom  part, the quantities
        ! have to be taken from their symmetric counterpart (2pi-theta).
        ! loop over all normal points of this rank
        do kd = 1,grid_out%loc_n_r
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
                    
                    ! set min
                    if (tb_sym) then
                        fund_theta_int_lo = -1
                        fund_theta_int_hi = 1
                    else
                        fund_theta_int_lo = 0
                        fund_theta_int_hi = 2
                    end if
                    
                    ! add or subtract  2pi to the poloidal angle until  it is in
                    ! the fundamental  range indicated by  fund_theta_int_lo and
                    ! hi
                    if (theta_loc.lt.fund_theta_int_lo*pi) then
                        do while (theta_loc.lt.fund_theta_int_lo*pi)
                            theta_loc = theta_loc + 2*pi
                            fund_theta_int_displ(id,jd,kd) = &
                                &fund_theta_int_displ(id,jd,kd) + 1
                        end do
                    else if (theta_loc.gt.fund_theta_int_hi*pi) then
                        do while (theta_loc.gt.fund_theta_int_hi*pi)
                            theta_loc = theta_loc - 2*pi
                            fund_theta_int_displ(id,jd,kd) = &
                                &fund_theta_int_displ(id,jd,kd) - 1
                        end do
                    end if
                end do
            end do
        end do
        
        ! clean up
        nullify(theta_in)
    end function get_ang_interp_data_HEL
    
    ! Interpolate  variables resulting  from HELENA  equilibria to  another grid
    ! (angularly).
    !   - equilibrium variables: flux variables (no need to convert) and derived
    !     quantities
    !   - metric variables: jac_FD, g_FD, h_FD
    !   - vectorial perturbation variables: U_i, DU_i
    !   - tensorial perturbation variables: PV_i, KV_i
    ! Also, a message can be printed if a grid name is passed.
    ! Note: By default the interpolated  quantities overwrite the original ones,
    ! but alternative output variables can be provided.
    ! Also  note   that  to  interpolate  the   metric  quantities,  equilibrium
    ! quantities  are needed  as well,  to set  up the  transformation matrices,
    ! through eq_met.
    integer function interp_HEL_on_grid(grid_eq,grid_eq_out,eq,met,X_1,X_2,&
        &eq_out,met_out,X_1_out,X_2_out,eq_met,grid_name) result(ierr)
        
        use num_vars, only: use_pol_flux_F
        use utilities, only: calc_mult, c
        
        character(*), parameter :: rout_name = 'interp_HEL_on_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq, grid_eq_out                     ! general and field-aligned equilibrium grid
        type(eq_type), intent(inout), optional :: eq                            ! general equilibrium variables
        type(met_type), intent(inout), optional :: met                          ! general metric variables
        type(X_1_type), intent(inout), optional  :: X_1                         ! general vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2                         ! general tensorial perturbation variables
        type(eq_type), intent(inout), optional :: eq_out                        ! field-aligned equilibrium variables
        type(met_type), intent(inout), optional :: met_out                      ! field-aligned metric variables
        type(X_1_type), intent(inout), optional  :: X_1_out                     ! field-aligned vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2_out                     ! field-aligned tensorial perturbation variables
        type(eq_type), intent(in), optional :: eq_met                           ! general equilibrium variables for metric interpolation
        character(len=*), intent(in), optional :: grid_name                     ! name of grid to which to adapt quantities
        
        ! local variables
        real(dp), allocatable :: theta_i(:,:,:)                                 ! interpolation index
        integer, allocatable :: fund_theta_int_displ(:,:,:)                     ! displacement of fundamental theta interval
        logical :: tb_sym                                                       ! whether there is top-bottom symmetry (relevant only for HELENA)
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: T_od(:)                                        ! nonzero off-diagonal part of T
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        if (present(grid_name)) then
            call writo('Adapting quantities to '//trim(grid_name))
            call lvl_ud(1)
        end if
        
        ! tests
        if (present(met).and..not.present(eq_met)) then
            ierr =1 
            err_msg = 'To interpolate metric variables, equilibrium has to be &
                &provided as well through eq_met'
            CHCKERR(err_msg)
        end if
        
        ! set up tb_sym
        if (ias.eq.0) then
            tb_sym = .true.
        else
            tb_sym = .false.
        end if
        
        ! get angular interpolation factors
        ierr = get_ang_interp_data_HEL(grid_eq,grid_eq_out,theta_i,&
            &fund_theta_int_displ,tb_sym=tb_sym,use_E=.false.)
        CHCKERR('')
        
        ! adapt equilibrium quantities
        if (present(eq)) then
            if (present(eq_out)) then
                eq_out%pres_FD = eq%pres_FD
                eq_out%q_saf_FD = eq%q_saf_FD
                eq_out%rot_t_FD = eq%rot_t_FD
                eq_out%flux_p_FD = eq%flux_p_FD
                eq_out%flux_t_FD = eq%flux_t_FD
                eq_out%rho = eq%rho
            !else
                ! do nothing
            end if
            if (present(eq_out)) then
                ierr = interp_var_3D_real(eq%S,theta_i,eq_out%S)
            else
                ierr = interp_var_3D_real_ow(eq%S,theta_i)
            end if
            CHCKERR('')
            if (present(eq_out)) then
                ierr = interp_var_3D_real(eq%sigma,theta_i,eq_out%sigma)
            else
                ierr = interp_var_3D_real_ow(eq%sigma,theta_i)
            end if
            CHCKERR('')
            if (present(eq_out)) then
                ierr = interp_var_3D_real(eq%kappa_n,theta_i,eq_out%kappa_n)
            else
                ierr = interp_var_3D_real_ow(eq%kappa_n,theta_i)
            end if
            CHCKERR('')
            if (present(eq_out)) then
                ierr = interp_var_3D_real(eq%kappa_g,theta_i,eq_out%kappa_g,&
                    &sym_type=2)
            else
                ierr = interp_var_3D_real_ow(eq%kappa_g,theta_i,sym_type=2)
            end if
            CHCKERR('')
        end if
        
        ! metric quantities
        if (present(met)) then
            ! setup T_od
            allocate(T_od(grid_eq%loc_n_r))
            if (use_pol_flux_F) then
                T_od = 2*pi*eq_met%q_saf_FD(:,0)
            else
                T_od = -2*pi*eq_met%rot_t_FD(:,0)
            end if
            ! interpolate
            if (present(met_out)) then
                ierr = interp_var_6D_real(met%jac_FD,theta_i,met_out%jac_FD)
            else
                ierr = interp_var_6D_real_ow(met%jac_FD,theta_i)
            end if
            CHCKERR('')
            if (present(met_out)) then
                ierr = interp_var_7D_real(met%g_FD,theta_i,met_out%g_FD,&
                    &sym_type=3,T_od=T_od,&
                    &fund_theta_int_displ=fund_theta_int_displ)
            else
                ierr = interp_var_7D_real_ow(met%g_FD,theta_i,sym_type=3,&
                    &T_od=T_od,fund_theta_int_displ=fund_theta_int_displ)
            end if
            CHCKERR('')
            if (present(met_out)) then
                ierr = interp_var_7D_real(met%h_FD,theta_i,met_out%h_FD,&
                    &sym_type=3,T_od=T_od,&
                    &fund_theta_int_displ=fund_theta_int_displ)
            else
                ierr = interp_var_7D_real_ow(met%h_FD,theta_i,sym_type=3,&
                    &T_od=T_od,fund_theta_int_displ=fund_theta_int_displ)
            end if
            CHCKERR('')
        end if
        
        ! adapt vectorial perturbation quantities
        if (present(X_1)) then
            if (present(X_1_out)) then
                ierr = interp_var_4D_complex(X_1%U_0,theta_i,X_1_out%U_0,&
                    &sym_type=2)
            else
                ierr = interp_var_4D_complex_ow(X_1%U_0,theta_i,sym_type=2)
            end if
            CHCKERR('')
            if (present(X_1_out)) then
                ierr = interp_var_4D_complex(X_1%U_1,theta_i,X_1_out%U_1,&
                    &sym_type=2)
            else
                ierr = interp_var_4D_complex_ow(X_1%U_1,theta_i,sym_type=2)
            end if
            CHCKERR('')
            if (present(X_1_out)) then
                ierr = interp_var_4D_complex(X_1%DU_0,theta_i,X_1_out%DU_0,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_1%DU_0,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_1_out)) then
                ierr = interp_var_4D_complex(X_1%DU_1,theta_i,X_1_out%DU_1,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_1%DU_1,theta_i,sym_type=1)
            end if
            CHCKERR('')
        end if
        
        ! adapt tensorial perturbation quantities
        if (present(X_2)) then
            if (present(X_2_out)) then
                X_2_out%vac_res = X_2%vac_res
            else
                ! do nothing
            end if
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%PV_0,theta_i,X_2_out%PV_0,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%PV_0,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%PV_1,theta_i,X_2_out%PV_1,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%PV_1,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%PV_2,theta_i,X_2_out%PV_2,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%PV_2,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%KV_0,theta_i,X_2_out%KV_0,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%KV_0,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%KV_1,theta_i,X_2_out%KV_1,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%KV_1,theta_i,sym_type=1)
            end if
            CHCKERR('')
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%KV_2,theta_i,X_2_out%KV_2,&
                    &sym_type=1)
            else
                ierr = interp_var_4D_complex_ow(X_2%KV_2,theta_i,sym_type=1)
            end if
            CHCKERR('')
        end if
        
        ! user output
        if (present(grid_name)) then
            call lvl_ud(-1)
            call writo('Quantities adapted to '//trim(grid_name))
        end if
    contains
        ! Interpolate a  variable defined  on an  axisymmetric grid  at poloidal
        ! angles indicated by the interpolation factors theta_i.
        ! There  is  an  optional  variable   sym_type  that  allows  for  extra
        ! operations to be done on the variable if top-down symmetry is applied:
        !   - sym_type = 0: var(-theta) = var(theta)
        !   - sym_type = 1: var(-theta) = var(theta)*
        !   - sym_type = 2: var(-theta) = -var(theta)*
        !   - sym_type = 3: var(-theta) = T var(theta) (T)^T
        ! where symmetry type 5 is only valid for tensorial quantities (7D). The
        ! transformation matrix, a function of  only the flux coordinate, has to
        ! be  passed as  well. This  is useful  for the  transformations of  the
        ! metric quantities as seen in [ADD REF].
        ! When top-down  symmetry has been  used to calculate  the interpolation
        ! factors, this is indicated by a  negative factor instead of a positive
        ! one.
        ! Note that this is also  correct  if i_hi = i_lo.
        ! Note also that for the overwrite  versions ('ow'), the lower limits of
        ! indices are 0  for the dimensions corresponding  to derivatives. These
        ! allocations are done here and persist into the calling functions.
        integer function interp_var_3D_real(var,theta_i,var_int,sym_type) &
            &result(ierr)                                                       ! 3D_real version, separate output variable
            
            ! input / output
            real(dp), intent(in) :: var(:,:,:)                                  ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:)                           ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd                                               ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! initialize ierr
            ierr = 0
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! test symmetry types
            if (sym_type_loc.lt.0 .or. sym_type_loc.gt.2) then
                ierr = 1
                err_msg = 'Nonsensible symmetry type '//&
                    &trim(i2str(sym_type_loc))
                CHCKERR(err_msg)
            end if
            
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
        end function interp_var_3D_real
        integer function interp_var_3D_real_ow(var,theta_i,sym_type) &
            &result(ierr)                                                       ! 3D_real version, overwriting variable
            
            ! input / output
            real(dp), intent(inout), allocatable :: var(:,:,:)                  ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            real(dp), allocatable :: var_3D_loc(:,:,:)                          ! local 3D variable
            integer :: var_dim(3,2)                                             ! dimensions of variable
            
            ! initialize ierr
            ierr = 0
            
            ! set up variable dimensions
            var_dim(:,1) = [1,1,1]
            var_dim(:,2) = [shape(theta_i)]
            
            ! set up local 3D variable
            allocate(var_3D_loc(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2)))
            
            ! call normal routine
            ierr = interp_var_3D_real(var,theta_i,var_3D_loc,sym_type)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2)))
            
            ! overwrite variable
            var = var_3D_loc
        end function interp_var_3D_real_ow
        integer function interp_var_4D_complex(var,theta_i,var_int,sym_type) &
            &result(ierr)                                                       ! 4D_complex version, separate output variable
            
            ! input / output
            complex(dp), intent(in) :: var(:,:,:,:)                             ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            complex(dp), intent(inout) :: var_int(:,:,:,:)                      ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd                                               ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! initialize ierr
            ierr = 0
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! test symmetry types
            if (sym_type_loc.lt.0 .or. sym_type_loc.gt.2) then
                ierr = 1
                err_msg = 'Nonsensible symmetry type '//&
                    &trim(i2str(sym_type_loc))
                CHCKERR(err_msg)
            end if
            
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
        end function interp_var_4D_complex
        integer function interp_var_4D_complex_ow(var,theta_i,sym_type) &
            &result(ierr)                                                       ! 4D_complex version, overwriting variable
            
            ! input / output
            complex(dp), intent(inout), allocatable :: var(:,:,:,:)             ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            complex(dp), allocatable :: var_4D_loc(:,:,:,:)                     ! local 4D variable
            integer :: var_dim(4,2)                                             ! dimensions of variable
            
            ! initialize ierr
            ierr = 0
            
            ! set up variable dimensions
            var_dim(:,1) = [1,1,1,1]
            var_dim(:,2) = [shape(theta_i),size(var,4)]
            
            ! set up local 4D variable
            allocate(var_4D_loc(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2)))
            
            ! call normal routine
            ierr = interp_var_4D_complex(var,theta_i,var_4D_loc,sym_type)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2)))
            
            ! overwrite variable
            var = var_4D_loc
        end function interp_var_4D_complex_ow
        integer function interp_var_6D_real(var,theta_i,var_int,sym_type) &
            &result(ierr)                                                       ! 6D_real version, separate output variable
            
            ! input / output
            real(dp), intent(in) :: var(:,:,:,:,:,:)                            ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:,:,:,:)                     ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd, ld                                           ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            
            ! initialize ierr
            ierr = 0
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! test symmetry types
            if (sym_type_loc.lt.0 .or. sym_type_loc.gt.2) then
                ierr = 1
                err_msg = 'Nonsensible symmetry type '//&
                    &trim(i2str(sym_type_loc))
                CHCKERR(err_msg)
            end if
            
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
                            ! adapt the derivatives
                            if (use_pol_flux_F) then
                                do ld = 1,size(var_int,6)
                                    var_int(id,jd,kd,:,:,ld) = &
                                        &(-1)**(ld+1)*var_int(id,jd,kd,:,:,ld)  ! index + 1 because not starting from 0 but from 1
                                end do
                            else
                                ierr = 1
                                err_msg = '!!! INTERP_VAR_6D_REAL NOT YET &
                                    &IMPLEMENTED FOR TOROIDAL FLUX !!!'
                                CHCKERR(err_msg)
                            end if
                        end if
                    end do
                end do
            end do
        end function interp_var_6D_real
        integer function interp_var_6D_real_ow(var,theta_i,sym_type) &
            &result(ierr)                                                       ! 6D_real version, overwriting variable
            
            ! input / output
            real(dp), intent(inout), pointer :: var(:,:,:,:,:,:)                ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            real(dp), allocatable :: var_6D_loc(:,:,:,:,:,:)                    ! local 6D variable
            integer :: var_dim(6,2)                                             ! dimensions of variable
            
            ! initialize ierr
            ierr = 0
            
            ! set up variable dimensions
            var_dim(:,1) = [1,1,1,0,0,0]
            var_dim(:,2) = [shape(theta_i),size(var,4)-1,size(var,5)-1,&
                &size(var,6)-1]
            
            ! set up local 6D variable
            allocate(var_6D_loc(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2),var_dim(5,1):var_dim(5,2),&
                &var_dim(6,1):var_dim(6,2)))
            
            ! call normal routine
            ierr = interp_var_6D_real(var,theta_i,var_6D_loc,sym_type)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2),var_dim(5,1):var_dim(5,2),&
                &var_dim(6,1):var_dim(6,2)))
            
            ! overwrite variable
            var = var_6D_loc
        end function interp_var_6D_real_ow
        integer function interp_var_7D_real(var,theta_i,var_int,sym_type,T_od,&
            &fund_theta_int_displ) result(ierr)                                 ! 7D_real version, separate output variable
            
            ! input / output
            real(dp), intent(in) :: var(:,:,:,:,:,:,:)                          ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            real(dp), intent(inout) :: var_int(:,:,:,:,:,:,:)                   ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            real(dp), intent(in), optional :: T_od(:)                           ! off-diagonal nonzero element of T in T var (T)^T
            integer, allocatable :: fund_theta_int_displ(:,:,:)                 ! displacement of fundamental theta interval
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd, ld                                           ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            real(dp), allocatable :: var_loc(:)                                 ! local variable
            integer :: n                                                        ! size of tensor matrix
            real(dp), allocatable :: T(:,:)                                     ! transformation matrix for metric factors
            
            ! initialize ierr
            ierr = 0
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! test symmetry types
            if (sym_type_loc.lt.0 .or. sym_type_loc.gt.3) then
                ierr = 1
                err_msg = 'Nonsensible symmetry type '//&
                    &trim(i2str(sym_type_loc))
                CHCKERR(err_msg)
            end if
            if (sym_type_loc.eq.3 .and. .not.present(T_od)) then
                ierr = 1
                err_msg = 'Need to pass nonzero off-diagonal element of &
                    &transformation matrix T for symmetry type 3'
                CHCKERR(err_msg)
            end if
            
            ! set up local variables for symmetry type 3
            if (sym_type_loc.eq.3) then
                allocate(T(size(T_od),9))                                       ! asymmetric matrix
                T = 0
                n = nint(sqrt(size(T,2)*1._dp))
                if (n**2.ne.size(T,2)) then
                    ierr = 1
                    err_msg = 'T does not have a proper size'
                    CHCKERR(err_msg)
                end if
                allocate(var_loc(n**2))
            end if
            
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
                        
                        ! apply symmetries
                        if (sym_type_loc.eq.2) then
                            if (theta_i(id,jd,kd).lt.0) &
                                &var_int(id,jd,kd,:,:,:,:) = &
                                &- var_int(id,jd,kd,:,:,:,:)
                        elseif (sym_type_loc.eq.3) then
                            ! apply possible transformation with T
                            if (fund_theta_int_displ(id,jd,kd).ne.0) then
                                T(:,c([2,1],.false.)) = &
                                    &T_od*fund_theta_int_displ(id,jd,kd)
                                ierr = calc_mult(T(kd,:),&
                                    &var_int(id,jd,kd,:,1,1,1),var_loc,n,&      ! index + 1 because not starting from 0 but 1
                                    &transp=[.false.,.false.])
                                CHCKERR('')
                                ierr = calc_mult(var_loc,T(kd,:),&
                                    &var_int(id,jd,kd,:,1,1,1),n,&
                                    &transp=[.false.,.true.])
                            end if
                        end if
                    end do
                end do
            end do
            
            if (use_pol_flux_F) then
                ! iterate over all alpha derivatives
                do ld = 2,size(var_int,5)                                       ! The indices start at 1, not at 0
                    var_int(:,:,:,:,ld,:,:) = 0._dp
                end do
                ! iterate over all psi derivatives
                write(*,*) '!!! NOT IMPLEMENTED !!!'
                ! iterate over all theta derivatives
                write(*,*) '!!! NOT IMPLEMENTED !!!'
            else
                ierr = 1
                err_msg = '!!! INTERP_VAR_7D_REAL NOT YET IMPLEMENTED FOR &
                    &TOROIDAL FLUX !!!'
                CHCKERR(err_msg)
            end if
            
            ! iterate over all normal points
            do kd = 1,size(var_int,3)
                ! iterate over all geodesical points
                do jd = 1,size(var_int,2)
                    ! iterate over all parallel points
                    do id = 1,size(var_int,1)
                        ! adapt the derivatives if TBS
                        if (theta_i(id,jd,kd).lt.0) then
                            if (use_pol_flux_F) then
                                do ld = 1,size(var_int,6)
                                    var_int(id,jd,kd,:,:,:,ld) = (-1)**(ld+1)*&
                                        &var_int(id,jd,kd,:,:,:,ld)             ! index + 1 because not starting from 0 but from 1
                                end do
                            else
                                ierr = 1
                                err_msg = '!!! INTERP_VAR_7D_REAL NOT YET &
                                    &IMPLEMENTED FOR TOROIDAL FLUX !!!'
                                CHCKERR(err_msg)
                            end if
                        end if
                    end do
                end do
            end do
        end function interp_var_7D_real
        integer function interp_var_7D_real_ow(var,theta_i,sym_type,T_od,&
            &fund_theta_int_displ) result(ierr)                                 ! 7D_real version, overwriting variable
            
            ! input / output
            real(dp), intent(inout), pointer :: var(:,:,:,:,:,:,:)              ! variable to be interpolated
            real(dp), intent(in) :: theta_i(:,:,:)                              ! angular coordinate theta at which to interpolate
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            real(dp), intent(in), optional :: T_od(:)                           ! off-diagonal nonzero element of T in T var (T)^T
            integer, allocatable :: fund_theta_int_displ(:,:,:)                 ! displacement of fundamental theta interval
            
            ! local variables
            real(dp), allocatable :: var_7D_loc(:,:,:,:,:,:,:)                  ! local 7D variable
            integer :: var_dim(7,2)                                             ! dimensions of variable
            
            ! initialize ierr
            ierr = 0
            
            ! set up variable dimensions
            var_dim(:,1) = [1,1,1,1,0,0,0]
            var_dim(:,2) = [shape(theta_i),size(var,4),size(var,5)-1,&
                &size(var,6)-1,size(var,7)-1]
            
            ! set up local 7D variable
            allocate(var_7D_loc(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2),var_dim(5,1):var_dim(5,2),&
                &var_dim(6,1):var_dim(6,2),var_dim(7,1):var_dim(7,2)))
            
            ! call normal routine
            ierr = interp_var_7D_real(var,theta_i,var_7D_loc,sym_type,T_od,&
                &fund_theta_int_displ)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2),var_dim(5,1):var_dim(5,2),&
                &var_dim(6,1):var_dim(6,2),var_dim(7,1):var_dim(7,2)))
            
            ! overwrite variable
            var = var_7D_loc
        end function interp_var_7D_real_ow
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
        use num_vars, only: rank, norm_disc_prec_eq
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
        real(dp), allocatable :: jac_alt(:,:)                                   ! alternative calculation for jac
        real(dp), allocatable :: h_H_11_alt(:,:,:)                              ! alternative calculation for upper metric factor 11
        real(dp), allocatable :: h_H_12_alt(:,:,:)                              ! alternative calculation for upper metric factor 12
        real(dp), allocatable :: h_H_33_alt(:,:,:)                              ! alternative calculation for upper metric factor 33
        character(len=max_str_ln) :: file_name                                  ! name of plot file
        character(len=max_str_ln) :: description                                ! description of plot
        integer :: r_min = 4                                                    ! first normal index that has meaning
        real(dp), allocatable :: tempvar(:,:,:,:)                               ! temporary variable
        
        ! initialize ierr
        ierr = 0
        
        if (rank.eq.0) then
            ! user output
            call writo('Checking consistency of metric factors')
            call lvl_ud(1)
            
            ! calculate  the  auxiliary quantities  Zchi,  zpsi,  Rchi and  Rpsi
            ! containing the derivatives as well as jac
            allocate(Rchi(nchi,n_r),Rpsi(nchi,n_r))
            allocate(Zchi(nchi,n_r),Zpsi(nchi,n_r))
            allocate(jac(nchi,n_r))
            
            do id = 1,nchi
                ierr = calc_deriv(R_H(id,:),Rpsi(id,:),flux_p_H/(2*pi),1,&
                    &norm_disc_prec_eq+1)
                CHCKERR('')
                ierr = calc_deriv(Z_H(id,:),Zpsi(id,:),flux_p_H/(2*pi),1,&
                    &norm_disc_prec_eq+1)
                CHCKERR('')
            end do
            do kd = 1,n_r
                ierr = calc_deriv(R_H(:,kd),Rchi(:,kd),chi_H,1,&
                    &norm_disc_prec_eq+1)
                CHCKERR('')
                ierr = calc_deriv(Z_H(:,kd),Zchi(:,kd),chi_H,1,&
                    &norm_disc_prec_eq+1)
                CHCKERR('')
            end do
            jac = Zchi*Rpsi - Zpsi*Rchi 
            
            ! calculate Jacobian differently
            allocate(jac_alt(nchi,n_r))
            do kd = 1,n_r
                jac_alt(:,kd) = qs(kd)/(h_H_33(:,kd)*RBphi(kd))
            end do
            
            ! output jacobian
            ! set some variables
            file_name = 'TEST_jac_H'
            description = 'Testing whether the HELENA Jacobian is consistent.'
            
            ! plot difference
            call plot_diff_HDF5(reshape(R_H(:,r_min:n_r)*jac(:,r_min:n_r),&
                &[nchi,1,n_r-r_min+1]),reshape(jac_alt(:,r_min:n_r),&
                &[nchi,1,n_r-r_min+1]),file_name,description=description,&
                &output_message=.true.)
            
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
                ierr = calc_deriv(RBphi,tempvar(id,1,:,1),flux_p_H/(2*pi),1,&
                    &norm_disc_prec_eq)
                CHCKERR('')
                ierr = calc_deriv(pres_H,tempvar(id,1,:,2),flux_p_H/(2*pi),1,&
                    &norm_disc_prec_eq)
                CHCKERR('')
                ierr = calc_deriv(qs/RBphi*h_H_11(id,:),tempvar(id,1,:,3),&
                    &flux_p_H/(2*pi),1,norm_disc_prec_eq)
                CHCKERR('')
            end do
            do kd = 1,n_r
                ierr = calc_deriv(qs(kd)/RBphi(kd)*h_H_12(:,kd),&
                    &tempvar(:,1,kd,4),chi_H,1,norm_disc_prec_eq)
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
