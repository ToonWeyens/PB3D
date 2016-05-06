!------------------------------------------------------------------------------!
!   Operations on HELENA variables                                             !
!------------------------------------------------------------------------------!
module HELENA_ops
#include <PB3D_macros.h>
    use num_vars, only: pi
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, max_str_ln
    use grid_vars, only: grid_type, disc_type, dealloc_disc
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type
    use HELENA_vars
    
    implicit none
    private
    public interp_HEL_on_grid
#if ldebug
    public test_metrics_H
#endif

contains
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
    ! The  displacement of  the  theta interval  towards  the fundamental  theta
    ! interval is also outputted. For asymmetric variables this is for example:
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
        use num_utilities, only: con2dis
        
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
    
    ! Interpolate  variables resulting  from HELENA  equilibria to  another grid
    ! (angularly).  The input  and  output grid  to be  provided  depend on  the
    ! quantities to be interpolated:
    !   - equilibrium variables: flux variables (no need to convert) and derived
    !     quantities (need equilibrium grid)
    !   - metric variables: jac_FD (need equilibrium grid)
    !   - vectorial perturbation variables: U_i, DU_i (need perturbation grid)
    !   - tensorial perturbation variables: PV_i, KV_i (need perturbation grid)
    ! Also, a message can be printed if a grid name is passed.
    ! Note: the  metric coefficients are  interpolated and then  compensated for
    ! the straight-field-line coordinates as in [ADD REFERENCE].
    ! Note: By default the interpolated  quantities overwrite the original ones,
    ! but alternative output variables can be provided.
    ! Note:  as  the  equilibrium  and   perturbation  grid  are  not  generally
    ! identical,  this routine  has to  be called  separately for  the variables
    ! tabulated in either grid.
    integer function interp_HEL_on_grid(grid_in,grid_out,eq_2,X_1,X_2,&
        &eq_2_out,X_1_out,X_2_out,eq_1,grid_name) result(ierr)
        
        use num_vars, only: use_pol_flux_F
        
        character(*), parameter :: rout_name = 'interp_HEL_on_grid'
        
        ! input / output
        type(grid_type), intent(in) :: grid_in, grid_out                        ! input and output grid
        type(eq_2_type), intent(inout), optional :: eq_2                        ! general metric equilibrium variables
        type(X_1_type), intent(inout), optional  :: X_1                         ! general vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2                         ! general tensorial perturbation variables
        type(eq_2_type), intent(inout), optional :: eq_2_out                    ! field-aligned metric equilibrium variables
        type(X_1_type), intent(inout), optional  :: X_1_out                     ! field-aligned vectorial perturbation variables
        type(X_2_type), intent(inout), optional  :: X_2_out                     ! field-aligned tensorial perturbation variables
        type(eq_1_type), intent(in), optional :: eq_1                           ! general flux equilibrium variables for metric interpolation
        character(len=*), intent(in), optional :: grid_name                     ! name of grid to which to adapt quantities
        
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
                X_2_out%vac_res = X_2%vac_res
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
        integer function interp_var_6D_real_ow(var,sym_type) result(ierr)       ! 6D_real version, overwriting variable
            ! input / output
            real(dp), intent(inout), pointer :: var(:,:,:,:,:,:)                ! variable to be interpolated
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
        integer function interp_var_7D_real_ow(var,sym_type) result(ierr)       ! 7D_real version, overwriting variable
            
            ! input / output
            real(dp), intent(inout), pointer :: var(:,:,:,:,:,:,:)              ! variable to be interpolated
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
        use grid_utilities, only: setup_deriv_data, apply_disc
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
        type(disc_type) :: norm_deriv_data                                      ! data for normal derivative
        type(disc_type) :: ang_deriv_data                                       ! data for angular derivative
        
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
            
            ! calculate derivatives
            ierr = setup_deriv_data(flux_p_H/(2*pi),norm_deriv_data,1,&
                &norm_disc_prec_eq+1)
            CHCKERR('')
            ierr = apply_disc(R_H,norm_deriv_data,Rpsi,2)
            CHCKERR('')
            ierr = apply_disc(Z_H,norm_deriv_data,Zpsi,2)
            CHCKERR('')
            ierr = setup_deriv_data(chi_H,ang_deriv_data,1,norm_disc_prec_eq+1)
            CHCKERR('')
            ierr = apply_disc(R_H,ang_deriv_data,Rchi,1)
            CHCKERR('')
            ierr = apply_disc(Z_H,ang_deriv_data,Zchi,1)
            CHCKERR('')
            
            ! calculate Jacobian
            jac = Zchi*Rpsi - Zpsi*Rchi 
            
            ! calculate Jacobian differently
            allocate(jac_alt(nchi,n_r))
            do kd = 1,n_r
                jac_alt(:,kd) = qs_H(kd)/(h_H_33(:,kd)*RBphi_H(kd))
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
            ierr = apply_disc(RBphi_H,norm_deriv_data,tempvar(1,1,:,1))
            CHCKERR('')
            ierr = apply_disc(pres_H,norm_deriv_data,tempvar(1,1,:,2))
            CHCKERR('')
            do id = 1,nchi
                tempvar(id,1,:,1:2) = tempvar(1,1,:,1:2)
                ierr = apply_disc(qs_H/RBphi_H*h_H_11(id,:),norm_deriv_data,&
                    &tempvar(id,1,:,3))
                CHCKERR('')
            end do
            do kd = 1,n_r
                ierr = apply_disc(qs_H(kd)/RBphi_H(kd)*h_H_12(:,kd),&
                    &ang_deriv_data,tempvar(:,1,kd,4))
                CHCKERR('')
            end do
            
            ! calculate pressure  balance in tempvar(1)
            !   mu_0 p' = F/(qR^2) (d/d1 (h_11 q/F) + d/d2 (h_12 q/F) + q F')
            do kd = 1,n_r
                tempvar(:,1,kd,1) = -RBphi_H(kd)*h_H_33(:,kd)/qs_H(kd) * &
                    &(tempvar(:,1,kd,1)*qs_H(kd) + tempvar(:,1,kd,3) + &
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
            
            ! clean up
            call dealloc_disc(norm_deriv_data)
            call dealloc_disc(ang_deriv_data)
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end if
    end function test_metrics_H
#endif
end module HELENA_ops
