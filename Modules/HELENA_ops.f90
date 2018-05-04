!------------------------------------------------------------------------------!
!> Operations on HELENA variables.
!------------------------------------------------------------------------------!
module HELENA_ops
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    use HELENA_vars
    
    implicit none
    private
    public read_HEL, interp_HEL_on_grid
#if ldebug
    public test_metrics_H, test_harm_cont_H
#endif

contains
    !> Reads the HELENA equilibrium data
    !!
    !! Adapted from HELENA routine \c IODSK.
    !! 
    !! The variales in the HELENA mapping file are globalized in two ways:
    !!  - X and  Y are normalized w.r.t. vacuum geometric  axis \c R_vac and
    !!  toroidal field at the geometric axis \c B_vac.
    !!  - <tt> R[m] = R_vac[m] (1 + eps X[]) </tt>,
    !!  - <tt> Z[m] = R_vac[m] eps Y[] </tt>.
    !! 
    !! The covariant  toroidal field  \c F_H,  \c pres_H  and poloidal  flux are
    !! normalized  w.r.t  magnetic axis  \c  R_m  and  total toroidal  field  at
    !! magnetic axis \c B_m:
    !!  - <tt> RBphi[Tm]     = F_H[] R_m[m] B_m[T] </tt>,
    !!  - <tt> pres[N/m^2]   = pres_H[] (B_m[T])^2/mu_0[N/A^2] </tt>,
    !!  - <tt> flux_p[Tm^2]  = 2pi (s[])^2 cpsurf[] B_m[T] (R_m[m])^2 </tt>. 
    !!
    !! The first  normalization type  is the  HELENA normalization,  whereas the
    !! second is the MISHKA normalization.
    !! 
    !! Everything is translated to MISHKA  normalization to make comparison with
    !! MISHKA simple. This is done using the factors:
    !!  - <tt> radius[] = a[m] / R_m[m] </tt>,
    !!  - <tt> eps[] = a[m] / R_vac[m] </tt>,
    !! so that the expressions become:
    !!  - <tt> R[m]          = radius[] (1/eps[] + X[])           R_m[m]</tt>,
    !!  - <tt> Z[m]          = radius[] Y[]                       R_m[m]</tt>,
    !!  - <tt> RBphi[Tm]     = F_H[]                     B_m[T]   R_m[m]</tt>,
    !!  - <tt> pres[N/m^2]   = pres_H[]                  B_m[T]^2 mu_0[N/A^2]^-1</tt>,
    !!  - <tt> flux_p[Tm^2]  = 2pi (s[])^2 cpsurf[]      B_m[T]   R_m[T]^2</tt>.
    !!
    !! Finally, in  HELENA, the total  current \c I,  the poloidal beta  and the
    !! density at the geometric axis can be prescribed through:
    !!  - <tt> XIAB = mu_0 I / (a_vac B_vac) </tt>,
    !!  - <tt> BETAP = 8 pi S \<p\> / (I^2 mu_0) </tt>,
    !!  - <tt> ZN </tt>,
    !!
    !! where <tt> a_vac = eps R_vac </tt> and \c B_vac are vacuum quantities,
    !! \c  S is  the 2-D  cross-sectional  area and  \<p\> is  the 2-D  averaged
    !! pressure.
    !!
    !! \note To translate this to the MISHKA normalization factors as
    !!  - <tt> R_m = (eps/radius) R_vac </tt>,
    !!  - <tt> B_m = B_vac / B0 </tt>,
    !!
    !! \note where
    !!  - \c radius is in the mapping file (12), as well as in the HELENA output
    !!  (20).
    !!  - \c  eps is in the  mapping file (12), as  well as in the  HELENA input
    !!  (10) and output (20).
    !!  - \c  B0 is in the  mapping file (12), as  well as in the  HELENA output
    !!  (20).
    !! \note Furthermore,  the density on axis  can be specified as  \c ZN0 from
    !! HELENA input (10). \n The other  variables should probably not be touched
    !! for consistency.
    !! \note Finally, the variables IAS and B0 are in the mapping file (12) only
    !! in patched versions. See \ref tutorial_eq.
    !!
    !! \return ierr
    integer function read_HEL(n_r_in,use_pol_flux_H) result(ierr)
        use num_vars, only: eq_name, eq_i, max_deriv, tol_zero
        use num_utilities, only: calc_int, spline
        use HELENA_vars, only: pres_H, q_saf_H, rot_t_H, flux_p_H, flux_t_H, &
            &nchi, chi_H, ias, RBphi_H, R_H, Z_H, RMtoG_H, BMtoG_H
        
        character(*), parameter :: rout_name = 'read_HEL'
        
        ! input / output
        integer, intent(inout) :: n_r_in                                        !< nr. of normal points in input grid
        logical, intent(inout) :: use_pol_flux_H                                !< .true. if HELENA equilibrium is based on pol. flux
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, kd                                                       ! counters
        integer :: nchi_loc                                                     ! local nchi
        integer :: max_loc_Z(2)                                                 ! location of maximum Z
        real(dp), allocatable :: s_H(:)                                         ! flux coordinate s
        real(dp), allocatable :: curj(:)                                        ! toroidal current
        real(dp), allocatable :: vx(:), vy(:)                                   ! R and Z of surface
        real(dp) :: diff_Dp                                                     ! Dp difference
        real(dp) :: tol_diff_Dp                                                 ! tolerance on diff_Dp
        real(dp) :: Dp0, Dpe                                                    ! derivative of pressure on axis and surface
        real(dp) :: Dq0, Dqe                                                    ! derivative of safety factor on axis and surface
        real(dp) :: Dj0, Dje                                                    ! derivative of toroidal current on axis and surface
        real(dp) :: DF0, DFe                                                    ! derivative of RBphi on axis and surface
        real(dp) :: radius, raxis                                               ! minor radius, major radius normalized w.r.t.  magnetic axis (i.e. R_m)
        real(dp) :: B0                                                          ! B_geo/B_m
        real(dp) :: eps                                                         ! inverse aspect ratio
        real(dp) :: cpsurf                                                      ! poloidal flux on surface
        real(dp) :: ellip                                                       ! ellipticity
        real(dp) :: tria                                                        ! triangularity
        real(dp) :: delta_RZ(2)                                                 ! R and Z range
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Reading data from HELENA output "' &
            &// trim(eq_name) // '"')
        call lvl_ud(1)
        
        ! set error messages
        err_msg = 'Failed to read HELENA output file. &
            &Does it output IAS and B0? See tutorial!'
        
        ! rewind input file
        rewind(eq_i)
        
        ! read mapped equilibrium from disk
        read(eq_i,*,IOSTAT=ierr) n_r_in,ias                                     ! nr. normal points (JS0), top-bottom symmetry
        CHCKERR(err_msg)
        n_r_in = n_r_in + 1                                                     ! HELENA outputs nr. of normal points - 1
        
        allocate(s_H(n_r_in))                                                   ! flux coordinate
        read(eq_i,*,IOSTAT=ierr) (s_H(kd),kd=1,n_r_in)                          ! it is squared below, after reading cpsurf
        CHCKERR(err_msg)
        
        allocate(q_saf_H(n_r_in,0:max_deriv+1))                                 ! safety factor
        read(eq_i,*,IOSTAT=ierr) (q_saf_H(kd,0),kd=1,n_r_in)
        CHCKERR(err_msg)
        read(eq_i,*,IOSTAT=ierr) dq0, Dqe                                       ! first point, last point
        CHCKERR(err_msg)
        ! Note: Dq0 an  Dqe are in the local coordinate  of the finite elements,
        ! and are therefore not used here: q_saf_H(:,1) will be overwritten
        q_saf_H(1,1) = Dq0
        read(eq_i,*,IOSTAT=ierr) (q_saf_H(kd,1),kd=2,n_r_in)                    ! second to last point (again)
        CHCKERR(err_msg)
        
        allocate(curj(n_r_in))                                                  ! toroidal current
        read(eq_i,*,IOSTAT=ierr) (curj(kd),kd=1,n_r_in)
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) Dj0,Dje                                        ! derivative of toroidal current at first and last point
        CHCKERR(err_msg)
        
        read(eq_i,*,IOSTAT=ierr) nchi                                           ! nr. poloidal points
        CHCKERR(err_msg)
        nchi_loc = nchi
        if (ias.ne.0) nchi = nchi + 1                                           ! extend grid to 2pi if asymmetric (so doubling a point!)
        
        allocate(chi_H(nchi))                                                   ! poloidal points
        read(eq_i,*,IOSTAT=ierr) (chi_H(id),id=1,nchi_loc)
        CHCKERR(err_msg)
        if (ias.ne.0) chi_H(nchi) = 2*pi
        
        allocate(h_H_11(nchi,n_r_in))                                           ! upper metric factor 1,1
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_11(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_in*nchi_loc)                                     ! (gem11)
        CHCKERR(err_msg)
        h_H_11(:,1) = tol_zero                                                  ! first normal point is not given, so set to almost zero
        if (ias.ne.0) h_H_11(nchi,:) = h_H_11(1,:)
        
        allocate(h_H_12(nchi,n_r_in))                                           ! upper metric factor 1,2
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_12(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_in*nchi_loc)                                     ! (gem12)
        CHCKERR(err_msg)
        do id = 1,nchi-1
            ierr = spline(s_H(2:n_r_in),h_H_12(id,2:n_r_in),s_H(1:1),&
                &h_H_12(id,1:1),deriv=0,extrap=.true.)
            CHCKERR('')
        end do
        if (ias.ne.0) h_H_12(nchi,:) = h_H_12(1,:)
        
        read(eq_i,*,IOSTAT=ierr) cpsurf, radius                                 ! poloidal flux on surface, minor radius
        CHCKERR(err_msg)
        
        allocate(h_H_33(nchi,n_r_in))                                           ! lower metric factor 3,3
        read(eq_i,*,IOSTAT=ierr) &
            &(h_H_33(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),id=&
            &nchi_loc+1,n_r_in*nchi_loc)                                        ! (gem33)
        
        read(eq_i,*,IOSTAT=ierr) raxis, B0                                      ! major radius
        CHCKERR(err_msg)
        h_H_33(:,1) = raxis**2                                                  ! first normal point is degenerate, major radius
        if (ias.ne.0) h_H_33(nchi,:) = h_H_33(1,:)
        h_H_33 = 1._dp/h_H_33                                                   ! inverse is given
        
        allocate(pres_H(n_r_in,0:max_deriv+1))                                  ! pressure profile
        read(eq_i,*,IOSTAT=ierr) (pres_H(kd,0),kd=1,n_r_in)
        CHCKERR(err_msg)
        read(eq_i,*,IOSTAT=ierr) Dp0, Dpe                                       ! first point, last point
        CHCKERR(err_msg)
        
        allocate(RBphi_H(n_r_in,0:max_deriv+1))                                 ! R B_phi (= F)
        read(eq_i,*,IOSTAT=ierr) (RBphi_H(kd,0),kd=1,n_r_in)
        CHCKERR(err_msg)
        ! Note: DF0 and DFe contain rubbish and are not used.
        read(eq_i,*,IOSTAT=ierr) DF0, DFe                                       ! derivatives of R B_phi on axis and surface
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
        
        allocate(R_H(nchi,n_r_in))                                              ! major radius R
        read(eq_i,*,IOSTAT=ierr) (R_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_in*nchi_loc)                                     ! (xout)
        CHCKERR(err_msg)
        
        allocate(Z_H(nchi,n_r_in))                                              ! height Z
        read(eq_i,*,IOSTAT=ierr) (Z_H(mod(id-1,nchi_loc)+1,(id-1)/nchi_loc+1),&
            &id=nchi_loc+1,n_r_in*nchi_loc)                                     ! (yout)
        CHCKERR(err_msg)
        
        ! transform to MISHKA normalization
        allocate(flux_p_H(n_r_in,0:max_deriv+1))  
        flux_p_H(:,0) = 2*pi*s_H**2 * cpsurf                                    ! rescale flux coordinate (HELENA uses psi_pol/2pi as flux)
        radius = radius * raxis                                                 ! global length normalization with R_m
        R_H = radius*(1._dp/eps + R_H)                                          ! local normalization with a
        Z_H = radius*Z_H                                                        ! local normalization with a
        R_H(:,1) = raxis                                                        ! first normal point is degenerate, major radius
        Z_H(:,1) = 0._dp                                                        ! first normal point is degenerate, 0
        if (ias.ne.0) R_H(nchi,:) = R_H(1,:)
        if (ias.ne.0) Z_H(nchi,:) = Z_H(1,:)
        
        ! set conversion factors
        RMtoG_H = radius/eps                                                    ! R_geo / R_mag
        BMtoG_H = B0                                                            ! B_geo / B_mag
        
        ! calculate ellipticity and triangularity
        max_loc_Z = maxloc(Z_H)
        delta_RZ(1) = maxval(R_H)-minval(R_H)
        if (ias.eq.0) then
            delta_RZ(2) = 2*maxval(Z_H)
        else
            delta_RZ(2) = maxval(Z_H)-minval(Z_H)
        end if
        ellip = delta_RZ(2)/delta_RZ(1)
        tria = ((maxval(R_H)+minval(R_H))*0.5_dp-&
            &R_H(max_loc_Z(1),max_loc_Z(2)))/&
            &(delta_RZ(2)*0.5_dp)
        call writo('Calculated ellipticity: '//trim(r2str(ellip)))
        call writo('Calculated triangularity: '//trim(r2str(tria)))
        
        ! calculate derivatives of fluxes
        allocate(flux_t_H(n_r_in,0:max_deriv+1))
        flux_p_H(:,1) = 2*pi
        flux_p_H(:,2:max_deriv+1) = 0._dp
        flux_t_H(:,1) = q_saf_H(:,0)*flux_p_H(:,1)
        ierr = calc_int(flux_t_H(:,1),flux_p_H(:,0)/(2*pi),flux_t_H(:,0))
        CHCKERR('')
        do kd = 2,max_deriv+1
            ierr = spline(flux_p_H(:,0)/(2*pi),flux_t_H(:,1),&
                &flux_p_H(:,0)/(2*pi),flux_t_H(:,kd),deriv=kd-1)
            CHCKERR('')
        end do
        
        ! calculate derivatives of other variables
        allocate(rot_t_H(n_r_in,0:max_deriv+1))                                 ! rotational transform
        rot_t_H(:,0) = 1._dp/q_saf_H(:,0)
        do kd = 1,max_deriv+1
            ierr = spline(flux_p_H(:,0)/(2*pi),pres_H(:,0),&
                &flux_p_H(:,0)/(2*pi),pres_H(:,kd),deriv=kd,bcs=[0,1],&
                &bcs_val=[0._dp,Dpe*pi/(flux_p_H(n_r_in,0)*s_H(n_r_in))])
            CHCKERR('')
            ierr = spline(flux_p_H(:,0)/(2*pi),RBphi_H(:,0),&
                &flux_p_H(:,0)/(2*pi),RBphi_H(:,kd),deriv=kd)
            CHCKERR('')
            ierr = spline(flux_p_H(:,0)/(2*pi),rot_t_H(:,0),&
                &flux_p_H(:,0)/(2*pi),rot_t_H(:,kd),deriv=kd)
            CHCKERR('')
            ierr = spline(flux_p_H(:,0)/(2*pi),q_saf_H(:,0),&
                &flux_p_H(:,0)/(2*pi),q_saf_H(:,kd),deriv=kd)
            CHCKERR('')
        end do
        
        ! check for consistency
        tol_diff_Dp = abs(pres_H(1,0)-pres_H(n_r_in,0))/n_r_in
        diff_Dp = s_H(1)*flux_p_H(n_r_in,0)/pi*pres_H(1,1)-Dp0
        if (abs(diff_Dp).gt.tol_diff_Dp) then
            call writo('Difference between pressure derivative on axis is &
                &large ('//trim(r2str(diff_Dp))//')')
            call writo('Maybe increase precision or number of normal points',&
                &warning=.true.)
        end if
        diff_Dp = s_H(n_r_in)*flux_p_H(n_r_in,0)/pi*pres_H(n_r_in,1)-Dpe
        if (abs(diff_Dp).gt.tol_diff_Dp) then
            call writo('Difference between pressure derivative at edge is &
                &large ('//trim(r2str(diff_Dp))//')')
            call writo('Maybe increase precision or number of normal points',&
                &warning=.true.)
        end if
        
        !!! To plot the cross-section
        !!call print_ex_2D(['cross_section'],'cross_section',Z_H,x=R_H,&
            !!&draw=.false.)
        !!call draw_ex(['cross_section'],'cross_section',n_r_in,1,.false.)
        !!call print_ex_2D(['cross_section_inv'],'cross_section_inv',&
            !!&transpose(Z_H),x=transpose(R_H),draw=.false.)
        !!call draw_ex(['cross_section_inv'],'cross_section_inv',nchi,1,.false.)
       
        ! HELENA always uses the poloidal flux
        use_pol_flux_H = .true.
        
        ! user output
        call writo('HELENA output given on '//trim(i2str(nchi))//&
            &' poloidal and '//trim(i2str(n_r_in))//' normal points')
        if (ias.eq.0) call writo('The equilibrium is top-bottom symmetric')
        call lvl_ud(-1)
        call writo('Data from HELENA output succesfully read')
    end function read_HEL
    
    !> Calculate interpolation factors for  angular interpolation in \c grid_out
    !! of quantities defined on \c grid_in.
    !! 
    !! This version is specific for  an input grid corresponding to axisymmetric
    !! variables with optional  top-down symmetry, as is the  case for variables
    !! resulting from HELENA equilibria.
    !!
    !! The  output  of  a 3-D  array  of  real  values  for the  poloidal  angle
    !! \f$\theta\f$ where the  floored integer of each value  indicates the base
    !! index of the interpolated value in the output grid and the modulus is the
    !! the fraction towards the next integer.
    !!
    !! The flag \c tb_sym indicates optionally that there is top-bottom symmetry
    !! as well as axisymmetry. When there is top-down symmetry, the variables in
    !! the lower  half (i.e. \f$-\pi  < \theta <  0\f$) are calculated  from the
    !! variables  in  the  upper  half  using the  symmetry  properties  of  the
    !! variables.  To indicate  this, the  sign of  the interpolation  factor is
    !! inverted to a negative value.
    !!
    !! The  displacement of  the theta  interval towards  the fundamental  theta
    !! interval is also outputted. For asymmetric variables this is for example:
    !!  -  1  for \f$-2\pi \ldots  0\f$
    !!  -  0  for \f$ 0    \ldots  2\pi\f$
    !!  - -1  for \f$ 2\pi \ldots  4\pi\f$
    !!
    !! etc.
    !!
    !! For symmetric variables, this is for example:
    !!  -  1  for \f$-3\pi \ldots -2\pi\f$,   with symmetry property
    !!  -  1  for \f$-2\pi \ldots -\pi\f$
    !!  -  0  for \f$-\pi  \ldots  0\f$,      with symmetry property
    !!  -  0  for \f$ 0    \ldots  \pi\f$
    !!  - -1  for \f$ \pi  \ldots  2\pi\f$,   with symmetry property
    !!  - -1  for \f$ 3\pi \ldots  3\pi\f$
    !! etc.
    !!
    !! By default,  the variables in the  Flux coord. system are  used, but this
    !! can be changed optionally with the flag \c "use_E.
    !!
    !! \return ierr
    integer function get_ang_interp_data_HEL(grid_in,grid_out,theta_i,&
        &fund_theta_int_displ,tb_sym,use_E) result(ierr)
        use num_utilities, only: con2dis
        
        character(*), parameter :: rout_name = 'get_ang_interp_data_HEL'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  !< input grid
        type(grid_type), intent(in) :: grid_out                                 !< output grid
        real(dp), allocatable, intent(inout) :: theta_i(:,:,:)                  !< interpolation index
        integer, allocatable, intent(inout) :: fund_theta_int_displ(:,:,:)      !< displacement of fundamental theta interval
        logical, intent(in), optional :: tb_sym                                 !< top-bottom symmetry
        logical, intent(in), optional :: use_E                                  !< whether E is used instead of F
        
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
            write(*,*) grid_in%loc_n_r, grid_out%loc_n_r
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
                theta_in => grid_in%theta_E(:,1,kd)                             ! axisymmetric grid should have only one toroidal point
            else
                theta_in => grid_in%theta_F(:,1,kd)                             ! axisymmetric grid should have only one toroidal point
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
                    ! fund_theta_int_hi
                    if (theta_loc.lt.fund_theta_int_lo*pi) then
                        do while (theta_loc.lt.fund_theta_int_lo*pi)
                            theta_loc = theta_loc + 2*pi
                            fund_theta_int_displ(id,jd,kd) = &
                                &fund_theta_int_displ(id,jd,kd) - 1             ! theta interval lies to the left of fund. theta interval
                        end do
                    else if (theta_loc.gt.fund_theta_int_hi*pi) then
                        do while (theta_loc.gt.fund_theta_int_hi*pi)
                            theta_loc = theta_loc - 2*pi
                            fund_theta_int_displ(id,jd,kd) = &
                                &fund_theta_int_displ(id,jd,kd) + 1             ! theta interval lies to the right of fund. theta interval
                        end do
                    end if
                    ierr = con2dis(abs(theta_loc),theta_i(id,jd,kd),theta_in)
                    CHCKERR('')
                    if (theta_loc.lt.0) theta_i(id,jd,kd) = - theta_i(id,jd,kd)
                end do
            end do
        end do
        
        ! clean up
        nullify(theta_in)
    end function get_ang_interp_data_HEL
    
    !> Interpolate  variables resulting from  HELENA equilibria to  another grid
    !! (angularly).
    !!
    !! The input  and  output grid  to be  provided  depend on  the
    !! quantities to be interpolated:
    !!   -  equilibrium  variables: flux  variables  (no  need to  convert)  and
    !!   derived quantities (need equilibrium grid)
    !!   - metric variables: jac_FD (need equilibrium grid)
    !!   - vectorial perturbation variables: U_i, DU_i (need perturbation grid)
    !!   - tensorial perturbation variables: PV_i, KV_i (need perturbation grid)
    !! 
    !! Also, a message can be printed if a grid name is passed.
    !! 
    !! \note
    !!   -# The  metric coefficients are  interpolated and then  compensated for
    !!   the straight-field-line coordinates as in \cite Weyens3D .
    !!   -# By default the interpolated  quantities overwrite the original ones,
    !!   but alternative output variables can be provided.
    !!   -#  as  the  equilibrium  and   perturbation  grid  are  not  generally
    !!   identical, this routine  has to be called separately  for the variables
    !!   tabulated in either grid.
    !!
    !! \return ierr
    integer function interp_HEL_on_grid(grid_in,grid_out,eq_2,X_1,X_2,&
        &eq_2_out,X_1_out,X_2_out,eq_1,grid_name) result(ierr)
        
        use num_vars, only: use_pol_flux_F
        
        character(*), parameter :: rout_name = 'interp_HEL_on_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in                                  !< input grid
        type(grid_type), intent(in) :: grid_out                                 !< output grid
        type(eq_2_type), intent(inout), optional :: eq_2                        !< general metric equilibrium variables
        type(X_1_type), intent(inout), optional  :: X_1                         !< general vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2                         !< general tensorial perturbation variables
        type(eq_2_type), intent(inout), optional :: eq_2_out                    !< field-aligned metric equilibrium variables
        type(X_1_type), intent(inout), optional  :: X_1_out                     !< field-aligned vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2_out                     !< field-aligned tensorial perturbation variables
        type(eq_1_type), intent(in), optional :: eq_1                           !< general flux equilibrium variables for metric interpolation
        character(len=*), intent(in), optional :: grid_name                     !< name of grid to which to adapt quantities
        
        ! local variables
        real(dp), allocatable :: theta_i(:,:,:)                                 ! interpolation index
        integer, allocatable :: fund_theta_int_displ(:,:,:)                     ! displacement of fundamental theta interval
        logical :: tb_sym                                                       ! whether there is top-bottom symmetry (relevant only for HELENA)
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        if (present(grid_name)) then
            call writo('Adapt quantities to '//trim(grid_name))
            call lvl_ud(1)
        end if
        
        ! tests
        if (present(eq_2).and..not.present(eq_1)) then
            ierr =1 
            err_msg = 'To interpolate metric equilibrium variables, flux &
                &equilibrium has to be provided as well through eq_1'
            CHCKERR(err_msg)
        end if
        if (present(eq_2) .and. (present(X_1).or.present(X_2))) then
            ierr = 1
            err_msg = 'Only the quantities on one type of grid can be &
                &interpolated in a single call'
            CHCKERR(err_msg)
        end if
        
        ! set up tb_sym
        if (ias.eq.0) then
            tb_sym = .true.
        else
            tb_sym = .false.
        end if
        
        ! get angular interpolation factors
        ierr = get_ang_interp_data_HEL(grid_in,grid_out,theta_i,&
            &fund_theta_int_displ,tb_sym=tb_sym,use_E=.false.)
        CHCKERR('')
        
        ! adapt equilibrium quantities
        if (present(eq_2)) then
            if (present(eq_2_out)) then
                ierr = interp_var_6D_real(eq_2%jac_FD,eq_2_out%jac_FD)
                CHCKERR('')
                ierr = interp_var_7D_real(eq_2%g_FD,eq_2_out%g_FD,sym_type=3)
                CHCKERR('')
                ierr = interp_var_7D_real(eq_2%h_FD,eq_2_out%h_FD,sym_type=4)
                CHCKERR('')
                ierr = interp_var_3D_real(eq_2%S,eq_2_out%S)
                CHCKERR('')
                ierr = interp_var_3D_real(eq_2%sigma,eq_2_out%sigma)
                CHCKERR('')
                ierr = interp_var_3D_real(eq_2%kappa_n,eq_2_out%kappa_n)
                CHCKERR('')
                ierr = interp_var_3D_real(eq_2%kappa_g,eq_2_out%kappa_g,&
                    &sym_type=2)
                CHCKERR('')
            else
                ierr = interp_var_6D_real_ow(eq_2%jac_FD)
                CHCKERR('')
                ierr = interp_var_7D_real_ow(eq_2%g_FD,sym_type=3)
                CHCKERR('')
                ierr = interp_var_7D_real_ow(eq_2%h_FD,sym_type=4)
                CHCKERR('')
                ierr = interp_var_3D_real_ow(eq_2%S)
                CHCKERR('')
                ierr = interp_var_3D_real_ow(eq_2%sigma)
                CHCKERR('')
                ierr = interp_var_3D_real_ow(eq_2%kappa_n)
                CHCKERR('')
                ierr = interp_var_3D_real_ow(eq_2%kappa_g,sym_type=2)
                CHCKERR('')
            end if
        end if
        
        ! adapt vectorial perturbation quantities
        if (present(X_1)) then
            if (present(X_1_out)) then
                ierr = interp_var_4D_complex(X_1%U_0,X_1_out%U_0,sym_type=2)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_1%U_1,X_1_out%U_1,sym_type=2)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_1%DU_0,X_1_out%DU_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_1%DU_1,X_1_out%DU_1,sym_type=1)
                CHCKERR('')
            else
                ierr = interp_var_4D_complex_ow(X_1%U_0,sym_type=2)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_1%U_1,sym_type=2)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_1%DU_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_1%DU_1,sym_type=1)
                CHCKERR('')
            end if
        end if
        
        ! adapt tensorial perturbation quantities
        if (present(X_2)) then
            if (present(X_2_out)) then
                ierr = interp_var_4D_complex(X_2%PV_0,X_2_out%PV_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_2%PV_1,X_2_out%PV_1,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_2%PV_2,X_2_out%PV_2,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_2%KV_0,X_2_out%KV_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_2%KV_1,X_2_out%KV_1,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex(X_2%KV_2,X_2_out%KV_2,sym_type=1)
                CHCKERR('')
            else
                ierr = interp_var_4D_complex_ow(X_2%PV_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_2%PV_1,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_2%PV_2,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_2%KV_0,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_2%KV_1,sym_type=1)
                CHCKERR('')
                ierr = interp_var_4D_complex_ow(X_2%KV_2,sym_type=1)
                CHCKERR('')
            end if
        end if
        
        ! user output
        if (present(grid_name)) then
            call lvl_ud(-1)
        end if
    contains
        ! Interpolate a  variable defined  on an  axisymmetric grid  at poloidal
        ! angles indicated by the interpolation factors theta_i.
        ! There  is  an  optional  variable   sym_type  that  allows  for  extra
        ! operations to be done on the variable if top-down symmetry is applied:
        !   - sym_type = 0: var(-theta) = var(theta)
        !   - sym_type = 1: var(-theta) = var(theta)*
        !   - sym_type = 2: var(-theta) = -var(theta)*
        !   - sym_type = 3: var(-theta) = +/-var(theta)* for g_FD (See note)
        !   - sym_type = 4: var(-theta) = +/-var(theta)* for h_FD (See note)
        ! where symmetry type 5 is only valid for tensorial quantities (7D). The
        ! transformation matrix, a function of  only the flux coordinate, has to
        ! be  passed as  well. This  is useful  for the  transformations of  the
        ! metric quantities as seen in [ADD REF].
        ! When top-down  symmetry has been  used to calculate  the interpolation
        ! factors, this is indicated by a  negative factor instead of a positive
        ! one.
        ! Note: This is also correct if i_hi = i_lo.
        ! Note:  For  the  overwrite  versions   ('ow'),  the  lower  limits  of
        ! indices are 0  for the dimensions corresponding  to derivatives. These
        ! allocations are done here and persist into the calling functions.
        ! Note: The  6D and 7D  overwrite versions take  a pointer instead  of a
        ! allocatable variable  because they  are to be  used for  jacobians and
        ! metric coefficients, which are defined as pointers.
        ! Note: symmetry type 3 for 7D variables (i.e. metric coefficients), the
        ! relation var(-theta)  = (-1)^(i+j) var(theta)*  where i and j  are the
        ! indices in g_ij  and h^ij. Furthermore, since  the metric coefficients
        ! are for  the Flux coordinate system,  they need to be  adapted for the
        ! secular coordinate theta or phi.
        !> \private
        integer function interp_var_3D_real(var,var_int,sym_type) result(ierr)  ! 3D_real version, separate output variable
            ! input / output
            real(dp), intent(in) :: var(1:,1:,1:)                               ! variable to be interpolated
            real(dp), intent(inout) :: var_int(1:,1:,1:)                        ! interpolated var
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
        !> \private
        integer function interp_var_3D_real_ow(var,sym_type) result(ierr)       ! 3D_real version, overwriting variable
            ! input / output
            real(dp), intent(inout), allocatable :: var(:,:,:)                  ! variable to be interpolated
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
            ierr = interp_var_3D_real(var,var_3D_loc,sym_type)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2)))
            
            ! overwrite variable
            var = var_3D_loc
        end function interp_var_3D_real_ow
        !> \private
        integer function interp_var_4D_complex(var,var_int,sym_type) &
            &result(ierr)                                                       ! 4D_complex version, separate output variable
            
            ! input / output
            complex(dp), intent(in) :: var(1:,1:,1:,1:)                         ! variable to be interpolated
            complex(dp), intent(inout) :: var_int(1:,1:,1:,1:)                  ! interpolated var
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
        !> \private
        integer function interp_var_4D_complex_ow(var,sym_type) result(ierr)    ! 4D_complex version, overwriting variable
            ! input / output
            complex(dp), intent(inout), allocatable :: var(:,:,:,:)             ! variable to be interpolated
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
            ierr = interp_var_4D_complex(var,var_4D_loc,sym_type)
            CHCKERR('')
            
            ! reallocate variable
            deallocate(var)
            allocate(var(var_dim(1,1):var_dim(1,2),&
                &var_dim(2,1):var_dim(2,2),var_dim(3,1):var_dim(3,2),&
                &var_dim(4,1):var_dim(4,2)))
            
            ! overwrite variable
            var = var_4D_loc
        end function interp_var_4D_complex_ow
        !> \private
        integer function interp_var_6D_real(var,var_int,sym_type) result(ierr)  ! 6D_real version, separate output variable
            ! input / output
            real(dp), intent(in) :: var(1:,1:,1:,0:,0:,0:)                      ! variable to be interpolated
            real(dp), intent(inout) :: var_int(1:,1:,1:,0:,0:,0:)               ! interpolated var
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
                                do ld = 0,size(var_int,6)-1
                                    var_int(id,jd,kd,:,:,ld) = &
                                        &(-1)**ld*var_int(id,jd,kd,:,:,ld)
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
        !> \private
        integer function interp_var_6D_real_ow(var,sym_type) result(ierr)       ! 6D_real version, overwriting variable
            ! input / output
            real(dp), intent(inout), allocatable :: var(:,:,:,:,:,:)            ! variable to be interpolated
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
            ierr = interp_var_6D_real(var,var_6D_loc,sym_type)
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
        !> \private
        integer function interp_var_7D_real(var,var_int,sym_type) result(ierr)  ! 7D_real version, separate output variable
            use num_utilities, only: calc_mult, c, add_arr_mult, derivs
            use num_vars, only: max_deriv
            
            ! input / output
            real(dp), intent(in) :: var(1:,1:,1:,1:,0:,0:,0:)                   ! variable to be interpolated
            real(dp), intent(inout) :: var_int(1:,1:,1:,1:,0:,0:,0:)            ! interpolated var
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
            ! local variables
            integer :: i_lo, i_hi                                               ! upper and lower index
            integer :: id, jd, kd, ld                                           ! counters
            integer :: sym_type_loc                                             ! local version of symmetry type
            real(dp), allocatable :: T_met(:,:,:,:,:,:)                         ! transformation matrix for metric elements
            integer, allocatable :: derivs_loc(:,:)                             ! array of all unique derivatives
            integer :: c_loc(2)                                                 ! local c
            
            ! initialize ierr
            ierr = 0
            
            ! set up local symmetry type
            sym_type_loc = 0
            if (present(sym_type)) sym_type_loc = sym_type
            
            ! test symmetry types
            if (sym_type_loc.lt.0 .or. sym_type_loc.gt.4) then
                ierr = 1
                err_msg = 'Nonsensible symmetry type '//&
                    &trim(i2str(sym_type_loc))
                CHCKERR(err_msg)
            end if
            
            ! interpolate
            do kd = 1,size(var_int,3)
                do jd = 1,size(var_int,2)
                    do id = 1,size(var_int,1)
                        i_lo = floor(abs(theta_i(id,jd,kd)))
                        i_hi = ceiling(abs(theta_i(id,jd,kd)))
                        
                        var_int(id,jd,kd,:,:,:,:) = var(i_lo,1,kd,:,:,:,:) + &
                            &(var(i_hi,1,kd,:,:,:,:)-var(i_lo,1,kd,:,:,:,:))*&
                            &(abs(theta_i(id,jd,kd))-i_lo)                      ! because i_hi - i_lo = 1
                    end do
                end do
            end do
            
            ! inversion of poloidal coordinate if TBS and symmetry type 2, 3, 4
            if (sym_type_loc.eq.2 .or. sym_type_loc.eq.3 .or. &
                &sym_type_loc.eq.4) then
                do kd = 1,size(var_int,3)
                    do jd = 1,size(var_int,2)
                        do id = 1,size(var_int,1)
                            if (theta_i(id,jd,kd).lt.0) then
                                ! invert value
                                if (sym_type_loc.eq.2) then
                                    var_int(id,jd,kd,:,:,:,:) = &
                                        &- var_int(id,jd,kd,:,:,:,:)
                                else if (sym_type_loc.eq.3 .or. &
                                    &sym_type_loc.eq.4) then                    ! metric coefficients: assuming 3D and symmetric
                                    c_loc(1) = c([2,1],.true.)
                                    c_loc(2) = c([3,2],.true.)
                                    var_int(id,jd,kd,c_loc(1),:,:,:) = &
                                        &- var_int(id,jd,kd,c_loc(1),:,:,:)
                                    var_int(id,jd,kd,c_loc(2),:,:,:) = &
                                        &- var_int(id,jd,kd,c_loc(2),:,:,:)
                                end if
                                ! invert derivatives
                                if (use_pol_flux_F) then                        ! invert parallel derivative
                                    do ld = 0,size(var_int,7)-1
                                        var_int(id,jd,kd,:,:,:,ld) = (-1)**ld*&
                                            &var_int(id,jd,kd,:,:,:,ld)
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
            end if
            
            ! set up displacement matrix if symmetry type 3 or 4
            if (sym_type_loc.eq.3 .or. sym_type_loc.eq.4) then
                allocate(T_met(size(var_int,1),size(var_int,2),&
                    &size(var_int,3),0:size(var_int,5)-1,&
                    &0:size(var_int,6)-1,0:size(var_int,7)-1))
                T_met = 0._dp
                if (use_pol_flux_F) then
                    do ld = 0,max_deriv
                        do kd = 1,size(theta_i,3)
                            T_met(:,:,kd,0,ld,0) = &
                                &2*pi*eq_1%q_saf_FD(kd,ld+1)
                        end do
                    end do
                else
                    do ld = 0,max_deriv
                        do kd = 1,size(theta_i,3)
                            T_met(:,:,kd,0,ld,0) = &
                                &-2*pi*eq_1%rot_t_FD(kd,ld+1)
                        end do
                    end do
                end if
                do kd = 1,size(theta_i,3)
                    do jd = 1,size(theta_i,2)
                        do id = 1,size(theta_i,1)
                            T_met(id,jd,kd,:,:,:) = T_met(id,jd,kd,:,:,:)*&
                                &fund_theta_int_displ(id,jd,kd)
                        end do
                    end do
                end do
            end if
            
            ! correct g_2i
            if (sym_type_loc.eq.3) then                                         ! lower metric elements h_FD
                if (use_pol_flux_F) then                                        ! poloidal flux coordinates
                    do id = 0,max_deriv
                        derivs_loc = derivs(id)
                        do jd = 1,size(derivs_loc,2)
                            ! g_22 -> g_22 + T_met * g_32
                            c_loc(1) = c([3,2],.true.)
                            c_loc(2) = c([2,2],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! g_32 -> g_32 + T_met * g_31
                            c_loc(1) = c([3,1],.true.)
                            c_loc(2) = c([3,2],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! g_22 -> g_22 + T_met * g_32
                            c_loc(1) = c([3,2],.true.)
                            c_loc(2) = c([2,2],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! g_12 -> g_12 + T_met * g_11
                            c_loc(1) = c([1,1],.true.)
                            c_loc(2) = c([1,2],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                        end do
                    end do
                else                                                            ! toroidal flux coordinates
                    ierr = 1
                    err_msg = '!!! INTERP_VAR_7D_REAL NOT YET &
                        &IMPLEMENTED FOR TOROIDAL FLUX !!!'
                    CHCKERR(err_msg)
                end if
            end if
            
            ! correct h_1j
            if (sym_type_loc.eq.4) then                                         ! upper metric elements h_FD
                if (use_pol_flux_F) then                                        ! poloidal flux coordinates
                    do id = 0,max_deriv
                        derivs_loc = derivs(id)
                        do jd = 1,size(derivs_loc,2)
                            ! h^11 -> h^11 - T_met * h^12
                            c_loc(1) = c([1,2],.true.)
                            c_loc(2) = c([1,1],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),-T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! h^12 -> h^12 - T_met * h^22
                            c_loc(1) = c([2,2],.true.)
                            c_loc(2) = c([1,2],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),-T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! h^11 -> h^11 - T_met * h^12
                            c_loc(1) = c([1,2],.true.)
                            c_loc(2) = c([1,1],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),-T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                            ! h^13 -> h^13 - T_met * h^23
                            c_loc(1) = c([2,3],.true.)
                            c_loc(2) = c([1,3],.true.)
                            ierr = add_arr_mult(&
                                &var_int(:,:,:,c_loc(1),:,:,:),-T_met,&
                                &var_int(:,:,:,c_loc(2),&
                                &derivs_loc(1,jd),derivs_loc(2,jd),&
                                &derivs_loc(3,jd)),derivs_loc(:,jd))
                            CHCKERR('')
                        end do
                    end do
                else                                                            ! toroidal flux coordinates
                    ierr = 1
                    err_msg = '!!! INTERP_VAR_7D_REAL NOT YET &
                        &IMPLEMENTED FOR TOROIDAL FLUX !!!'
                    CHCKERR(err_msg)
                end if
            end if
        end function interp_var_7D_real
        !> \private
        integer function interp_var_7D_real_ow(var,sym_type) result(ierr)       ! 7D_real version, overwriting variable
            
            ! input / output
            real(dp), intent(inout), allocatable :: var(:,:,:,:,:,:,:)          ! variable to be interpolated
            integer, intent(in), optional :: sym_type                           ! optionally another type of symmetry
            
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
            ierr = interp_var_7D_real(var,var_7D_loc,sym_type)
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
    !> Checks whether the metric elements provided by HELENA are consistent with
    !! a direct calculation using the coordinate transformations \cite Weyens3D.
    !!
    !! Direct calculations used:
    !! \f[\begin{aligned}
    !!   \left|\nabla \psi\right|^2             &= \frac{1}{\mathcal{J}^2}
    !!      \left(\left(\frac{\partial Z}{\partial \chi}\right)^2 +
    !!      \left(\frac{\partial R}{\partial \chi}\right)^2\right) \\
    !!   \left|\nabla \psi \cdot \nabla \chi\right|   &= \frac{1}{\mathcal{J}^2}
    !!      \left(
    !!      \frac{\partial Z}{\partial \chi} \frac{\partial Z}{\partial \psi} - 
    !!      \frac{\partial R}{\partial \chi} \frac{\partial R}{\partial \psi}
    !!      \right) \\
    !!   \left|\nabla \chi\right|^2             &= \frac{1}{\mathcal{J}^2}
    !!      \left(\left(\frac{\partial Z}{\partial \psi}\right)^2 +
    !!      \left(\frac{\partial R}{\partial \psi}\right)^2\right) \\
    !!   \left|\nabla \phi\right|^2             &= \frac{1}{R^2}
    !! \end{aligned}\f]
    !! with
    !! \f[\mathcal{J} = 
    !!      \frac{\partial Z}{\partial \psi} \frac{\partial R}{\partial \chi} - 
    !!      \frac{\partial R}{\partial \psi} \frac{\partial Z}{\partial \chi}\f]
    !! 
    !! Also, test whether the pressure balance
    !! \f$\nabla p = \vec{J}\times\vec{B} \f$ is satisfied.
    !!
    !! \ldebug
    !!
    !! \return ierr
    integer function test_metrics_H() result(ierr)
        use num_vars, only: rank, norm_disc_prec_eq
        use output_ops, only: plot_diff_HDF5
        use input_utilities, only: get_int
        use grid_vars, only: n_r_eq
        use num_utilities, only: spline
        
        character(*), parameter :: rout_name = 'test_metrics_H'
        
        ! local variables
        integer :: id, kd                                                       ! counters
        integer :: bcs(2,3)                                                     ! boundary conditions (theta(even), theta(odd), r)
        real(dp) :: bcs_val(2,3)                                                ! values for boundary conditions
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
            allocate(Rchi(nchi,n_r_eq),Rpsi(nchi,n_r_eq))
            allocate(Zchi(nchi,n_r_eq),Zpsi(nchi,n_r_eq))
            allocate(jac(nchi,n_r_eq))
            
            ! set up boundary conditions
            if (ias.eq.0) then                                                  ! top-bottom symmetric
                bcs(:,1) = [1,1]                                                ! theta(even): zero first derivative
                bcs(:,2) = [2,2]                                                ! theta(odd): zero first derivative
            else
                bcs(:,1) = [-1,-1]                                              ! theta(even): periodic
                bcs(:,2) = [-1,-1]                                              ! theta(odd): periodic
            end if
            bcs(:,3) = [0,0]                                                    ! r: not-a-knot
            bcs_val = 0._dp
            
            ! calculate derivatives
            do id = 1,nchi
                ierr = spline(flux_p_H(:,0)/(2*pi),R_H(id,:),&
                    &flux_p_H(:,0)/(2*pi),Rpsi(id,:),ord=norm_disc_prec_eq,&
                    &deriv=1,bcs=bcs(:,3),bcs_val=bcs_val(:,3))
                CHCKERR('')
                ierr = spline(flux_p_H(:,0)/(2*pi),Z_H(id,:),&
                    &flux_p_H(:,0)/(2*pi),Zpsi(id,:),ord=norm_disc_prec_eq,&
                    &deriv=1,bcs=bcs(:,3),bcs_val=bcs_val(:,3))
                CHCKERR('')
            end do
            do kd = 1,n_r_eq
                ierr = spline(chi_H,R_H(:,kd),chi_H,Rchi(:,kd),&
                    &ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,1),&
                    &bcs_val=bcs_val(:,1))                                      ! even
                CHCKERR('')
                ierr = spline(chi_H,Z_H(:,kd),chi_H,Zchi(:,kd),&
                    &ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,2),&
                    &bcs_val=bcs_val(:,2))                                      ! odd
                CHCKERR('')
            end do
            
            ! calculate Jacobian
            do kd = 1,n_r_eq
                jac(:,kd) = q_saf_H(kd,0)/(h_H_33(:,kd)*RBphi_H(kd,0))
            end do
            
            ! calculate Jacobian differently
            allocate(jac_alt(nchi,n_r_eq))
            jac_alt = R_H*(Zchi*Rpsi - Zpsi*Rchi)
            
            ! output jacobian
            ! set some variables
            file_name = 'TEST_jac_H'
            description = 'Testing whether the HELENA Jacobian is consistent.'
            
            ! plot difference
            call plot_diff_HDF5(&
                &reshape(jac(:,r_min:n_r_eq),&
                &[nchi,1,n_r_eq-r_min+1]),reshape(jac_alt(:,r_min:n_r_eq),&
                &[nchi,1,n_r_eq-r_min+1]),file_name,descr=description,&
                &output_message=.true.)
            
            ! calculate the metric factors directly
            allocate(h_H_11_alt(nchi,1,n_r_eq))
            allocate(h_H_12_alt(nchi,1,n_r_eq))
            allocate(h_H_33_alt(nchi,1,n_r_eq))
            
            h_H_11_alt(:,1,:) = (R_H/jac)**2 * (Zchi**2 + Rchi**2)
            h_H_12_alt(:,1,:) = -(R_H/jac)**2 * (Zchi*Zpsi + Rchi*Rpsi)
            h_H_33_alt(:,1,:) = 1._dp/(R_H**2)
            
            ! output h_H_11
            ! set some variables
            file_name = 'TEST_h_H_1_1'
            description = 'Testing whether the HELENA metric factor h_1_1 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_11_alt(:,:,r_min:n_r_eq),&
                &reshape(h_H_11(:,r_min:n_r_eq),[nchi,1,n_r_eq-r_min+1]),&
                &file_name,descr=description,output_message=.true.)
            
            ! output h_H_12
            ! set some variables
            file_name = 'TEST_h_H_1_2'
            description = 'Testing whether the HELENA metric factor h_1_2 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_12_alt(:,:,r_min:n_r_eq),&
                &reshape(h_H_12(:,r_min:n_r_eq),[nchi,1,n_r_eq-r_min+1]),&
                &file_name,descr=description,output_message=.true.)
            
            ! output h_H_33
            ! set some variables
            file_name = 'TEST_h_H_3_3'
            description = 'Testing whether the HELENA metric factor h_3_3 is &
                &consistent.'
            
            ! plot difference
            call plot_diff_HDF5(h_H_33_alt(:,:,r_min:n_r_eq),&
                &reshape(h_H_33(:,r_min:n_r_eq),[nchi,1,n_r_eq-r_min+1]),&
                &file_name,descr=description,output_message=.true.)
            
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
            allocate(tempvar(nchi,1,n_r_eq,4))
            do kd = 1,n_r_eq
                tempvar(:,1,kd,1) = RBphi_H(kd,1)
                tempvar(:,1,kd,2) = pres_H(kd,1)
            end do
            do id = 1,nchi
                ierr = spline(flux_p_H(:,0)/(2*pi),&
                    &q_saf_H(:,0)/RBphi_H(:,0)*h_H_11(id,:),&
                    &flux_p_H(:,0)/(2*pi),tempvar(id,1,:,3),&
                    &ord=norm_disc_prec_eq,deriv=1,bcs=bcs(:,3),&
                    &bcs_val=bcs_val(:,3))
                CHCKERR('')
            end do
            do kd = 1,n_r_eq
                ierr = spline(chi_H,q_saf_H(kd,0)/RBphi_H(kd,0)*h_H_12(:,kd),&
                    &chi_H,tempvar(:,1,kd,4),ord=norm_disc_prec_eq,deriv=1,&
                    &bcs=bcs(:,2),bcs_val=bcs_val(:,2))                         ! odd
                CHCKERR('')
            end do
            
            ! calculate pressure  balance in tempvar(1)
            !   mu_0 p' = F/(qR^2) (d/d1 (h_11 q/F) + d/d2 (h_12 q/F) + q F')
            do kd = 1,n_r_eq
                tempvar(:,1,kd,1) = &
                    &-RBphi_H(kd,0)*h_H_33(:,kd)/q_saf_H(kd,0) * &
                    &(tempvar(:,1,kd,1)*q_saf_H(kd,0) + tempvar(:,1,kd,3) + &
                    &tempvar(:,1,kd,4))
            end do
            
            ! output difference with p'
            ! set some variables
            file_name = 'TEST_p_H'
            description = 'Testing whether the HELENA pressure balance is &
                &consistent'
            
            ! plot difference
            call plot_diff_HDF5(tempvar(:,:,:,1),tempvar(:,:,:,2),&
                &file_name,descr=description,output_message=.true.)
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end if
    end function test_metrics_H
    
    !> Investaige harmonic content of the HELENA variables.
    !!
    !! \note Run with one process.
    integer function test_harm_cont_H() result(ierr)
        use grid_vars, only: n_r_eq
        use grid_utilities, only: nufft
        use HELENA_vars, only: nchi, chi_H, ias, R_H, Z_H
        
        character(*), parameter :: rout_name = 'test_metrics_H'
        
        ! local variables
        integer :: kd                                                           ! counter
        integer :: nchi_full                                                    ! nr. of pol. points in full grid, without doubling
        real(dp), allocatable :: f_loc(:,:)                                     ! local Fourier modes
        real(dp), allocatable :: chi_full(:)                                    ! chi_H in full poloidal grid
        real(dp), allocatable :: R_H_full(:,:)                                  ! R_H in full poloidal grid
        real(dp), allocatable :: Z_H_full(:,:)                                  ! Z_H in full poloidal grid
        real(dp), allocatable :: R_H_F(:,:,:,:)                                 ! Fourier modes of R_H
        real(dp), allocatable :: Z_H_F(:,:,:,:)                                 ! Fourier modes of Z_H
        character(len=max_str_ln) :: plot_name                                  ! name of plot
        character(len=2) :: loc_str(2)                                          ! local string
        
        ! initialize ierr
        ierr = 0
        
        ! set up full R and Z
        if (ias.eq.0) then
            nchi_full = 2*(nchi-1)
            allocate(chi_full(nchi_full))
            allocate(R_H_full(nchi_full,n_r_eq))
            allocate(Z_H_full(nchi_full,n_r_eq))
            chi_full(1:nchi-1) =              chi_H(1:nchi-1)
            chi_full(2*(nchi-1):nchi:-1) =    2*pi - chi_H(2:nchi)
            R_H_full(1:nchi-1,:) =            R_H(1:nchi-1,:)
            R_H_full(2*(nchi-1):nchi:-1,:) =  R_H(2:nchi,:)
            Z_H_full(1:nchi-1,:) =            Z_H(1:nchi-1,:)
            Z_H_full(2*(nchi-1):nchi:-1,:) = -Z_H(2:nchi,:)
        else
            nchi_full = nchi-1
            allocate(chi_full(nchi_full))
            allocate(R_H_full(nchi_full,n_r_eq))
            allocate(Z_H_full(nchi_full,n_r_eq))
            chi_full(1:nchi_full)   = chi_H(1:nchi_full)
            R_H_full(1:nchi_full,:) = R_H(1:nchi_full,:)
            Z_H_full(1:nchi_full,:) = Z_H(1:nchi_full,:)
        end if
        
        ! calculate NUFFT
        do kd = 1,n_r_eq
            ierr = nufft(chi_full,R_H_full(:,kd),f_loc)
            CHCKERR('')
            if (kd.eq.1) allocate(R_H_F(size(f_loc,1),1,n_r_eq,2))
            R_H_F(:,1,kd,:) = f_loc
            ierr = nufft(chi_full,Z_H_full(:,kd),f_loc)
            CHCKERR('')
            if (kd.eq.1) allocate(Z_H_F(size(f_loc,1),1,n_r_eq,2))
            Z_H_F(:,1,kd,:) = f_loc
        end do
        
        ! plot
        loc_str(1) = 'Re'
        loc_str(2) = 'Im'
        do kd = 1,2
            plot_name = loc_str(kd)//'_R_H_F'
            call print_ex_2D(['1'],plot_name,R_H_F(:,1,:,kd),draw=.false.)
            call draw_ex(['1'],plot_name,n_r_eq,1,.false.)
            plot_name = loc_str(kd)//'_Z_H_F'
            call print_ex_2D(['1'],plot_name,Z_H_F(:,1,:,kd),draw=.false.)
            call draw_ex(['1'],plot_name,n_r_eq,1,.false.)
            
            plot_name = loc_str(kd)//'_R_H_F_log'
            call print_ex_2D(['1'],plot_name,&
                &log10(max(1.E-10_dp,abs(R_H_F(:,1,:,kd)))),&
                &draw=.false.)
            call draw_ex(['1'],plot_name,n_r_eq,1,.false.)
            plot_name = loc_str(kd)//'_Z_H_F_log'
            call print_ex_2D(['1'],plot_name,&
                &log10(max(1.E-10_dp,abs(Z_H_F(:,1,:,kd)))),&
                &draw=.false.)
            call draw_ex(['1'],plot_name,n_r_eq,1,.false.)
        end do
        call plot_HDF5(['real','imag'],'R_H_F',R_H_F)
        call plot_HDF5(['real','imag'],'R_H_F_log',&
            &log10(max(1.E-10_dp,abs(R_H_F))))
        call plot_HDF5(['real','imag'],'Z_H_F',Z_H_F)
        call plot_HDF5(['real','imag'],'Z_H_F_log',&
            &log10(max(1.E-10_dp,abs(Z_H_F))))
    end function
#endif
end module HELENA_ops
