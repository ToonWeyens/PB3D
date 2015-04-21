!------------------------------------------------------------------------------!
!   Operations  on the  solution  vectors  such as  returning  X,  U or  their !
!   derivates, etc.                                                            !
!------------------------------------------------------------------------------!
module sol_ops
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, pi
    use grid_vars, only: grid_type
    use eq_vars, only: eq_type
    use metric_vars, only: metric_type
    use X_vars, only: X_type

    implicit none
    private
    public calc_real_XUQ, plot_X_vecs
    
    ! interfaces
    interface calc_real_XUQ
        module procedure calc_real_XUQ_arr, calc_real_XUQ_ind
    end interface

contains
    ! calculates the  normal or geodesic  component, or optionally  the parallel
    ! derivative of the plasma perturbation  or the normal or geodesic component
    ! of the magnetic field perturbation in the perturbation grid for a range in
    ! normalized  time (1  corresponds to  one period).  The variable  XUQ_style
    ! determines which one:
    !   - XUQ_style = 1: X (supports parallel derivative)
    !   - XUQ_style = 2: U (supports parallel derivative)
    !   - XUQ_style = 3: Qn
    !   - XUQ_style = 4: Qg
    ! For Qn  and Qg,  the metric  variables have  to be  provided as  well. The
    ! output is given in the X grid for the normal part. The angular part of the
    ! output is given by  the equilibrium grid. These grids need  to be given in
    ! Flux coordinates.
    ! Note: The angular part of the X grid is neglected as it is assumed to be a
    ! 1D grid.
    integer function calc_real_XUQ_arr(grid_eq,eq,grid_X,X,X_id,XUQ_style,time,&
        &XUQ,met,deriv) result(ierr)                                            ! (time) array version
        use num_vars, only: use_pol_flux_F, grp_rank, grp_n_procs
        use utilities, only: con2dis, calc_deriv
        
        character(*), parameter :: rout_name = 'calc_real_XUQ_arr'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time(:)                                         ! time range
        real(dp), intent(inout) :: XUQ(:,:,:,:)                                 ! X, U, Qn or Qg
        type(metric_type), optional, intent(in) :: met                          ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        integer :: n_t                                                          ! number of time points
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id, jd, kd                                                   ! counter
        complex(dp) :: sqrt_X_val_norm                                          ! normalized sqrt(X_val)
        logical :: deriv_loc                                                    ! local copy of deriv
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        complex(dp), allocatable :: DX_vec(:)                                   ! normal derivative of X_vec for a specific mode
        complex(dp), allocatable :: fac_0(:,:,:), fac_1(:,:,:)                  ! factor to be multiplied with X and DX
        complex(dp), allocatable :: par_fac(:)                                  ! multiplicative factor due to parallel derivative
        complex(dp), allocatable :: XUQ_loc(:,:)                                ! complex XUQ without time at a normal point
        real(dp), allocatable :: grp_r_eq(:)                                    ! grp_r_F of X grid interpolated in eq grid
        integer :: i_lo, i_hi                                                   ! upper and lower index for interpolation of eq grid to X grid
        
        ! initialize ierr
        ierr = 0
        
        ! set up n_t
        n_t = size(time)
        
        ! tests
        if (grid_eq%n(1).ne.size(XUQ,1) .or. grid_eq%n(2).ne.size(XUQ,2) .or. &
            &grid_X%grp_n_r.ne.size(XUQ,3) .or. n_t.ne.size(XUQ,4)) then
            ierr = 1
            err_msg = 'XUQ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        if ((XUQ_style.eq.3 .or. XUQ_style.eq.4).and..not.present(met)) then
            ierr = 1
            err_msg = 'When calculating Q, also metric variables have to be &
                &provided'
            CHCKERR(err_msg)
        end if
        
        ! set up local grp_n_r of X grid
        if (grp_rank.lt.grp_n_procs-1) then
            grp_n_r_X_loc = grid_X%grp_n_r - 1                                  ! grp_n_r of X grid has ghost region
        else
            grp_n_r_X_loc = grid_X%grp_n_r
        end if
        
        ! set up local copy of deriv
        deriv_loc = .false.
        if (present(deriv)) deriv_loc = deriv
        
        ! set up normalized sqrt(X_val)
        sqrt_X_val_norm = sqrt(X%val(X_id))
        if (imagpart(sqrt_X_val_norm).gt.0) &
            &sqrt_X_val_norm = - sqrt_X_val_norm                                ! exploding solution, not the decaying one
        sqrt_X_val_norm = sqrt_X_val_norm / abs(sqrt_X_val_norm)                ! normalize
        
        ! calculate  normal  interpolation  tables for  equilibrium  grid  (this
        ! concerns eq variables and met variables)
        allocate(grp_r_eq(grid_X%grp_n_r))
        do kd = 1,grid_X%grp_n_r
            ierr = con2dis(grid_X%grp_r_F(kd),grp_r_eq(kd),grid_eq%grp_r_F)
            CHCKERR('')
        end do
        
        ! initialize multiplicative factors fac_0 and fac_0 and par_fac
        allocate(fac_0(grid_eq%n(1),grid_eq%n(2),grid_eq%grp_n_r))
        allocate(fac_1(grid_eq%n(1),grid_eq%n(2),grid_eq%grp_n_r))
        allocate(par_fac(grid_eq%grp_n_r))
        
        ! set up DX_vec
        allocate(DX_vec(grp_n_r_X_loc))
        
        ! initialize XUQ
        XUQ = 0._dp
        
        ! initialize XUQ_loc
        allocate(XUQ_loc(size(XUQ,1),size(XUQ,2)))
        
        ! iterate over all modes
        do jd = 1,X%n_mod
            ! set up parallel multiplicative factor par_fac in normal eq grid
            if (deriv_loc) then                                                 ! parallel derivative
                if (use_pol_flux_F) then
                    par_fac = iu*(X%n(jd)*eq%q_saf_FD(:,0)-X%m(jd))
                else
                    par_fac = iu*(X%n(jd)-X%m(jd)*eq%rot_t_FD(:,0))
                end if
            else
                par_fac = 1._dp
            end if
            
            ! set up multiplicative factors in eq grid
            select case (XUQ_style)
                case (1)                                                        ! calculating X
                    if (deriv_loc) then                                         ! parallel derivative
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = par_fac(kd)
                        end do
                    else
                        fac_0 = 1._dp
                    end if
                    fac_1 = 0._dp
                case (2)                                                        ! calculating U
                    if (deriv_loc) then                                         ! parallel derivative
                        fac_0 = X%DU_0(:,:,:,jd)
                        fac_1 = X%DU_1(:,:,:,jd)
                    else
                        fac_0 = X%U_0(:,:,:,jd)
                        fac_1 = X%U_1(:,:,:,jd)
                    end if
                case (3)                                                        ! calculating Qn
                    if (deriv_loc) then                                         ! parallel derivative
                        ierr = 1
                        err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                            &' calculation of parallel derivative not supported'
                        CHCKERR('')
                    else
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = par_fac(kd)/met%jac_FD(:,:,kd,0,0,0)
                        end do
                    end if
                    fac_1 = 0._dp
                case (4)                                                        ! calculating Qg
                    if (deriv_loc) then                                         ! parallel derivative
                        ierr = 1
                        err_msg = 'For XUQ_style = '//trim(i2str(XUQ_style))//&
                            &' calculation of parallel derivative not supported'
                        CHCKERR('')
                    else
                        do kd = 1,grid_eq%grp_n_r
                            fac_0(:,:,kd) = (X%DU_0(:,:,kd,jd)-&
                                &X%extra1(:,:,kd))/met%jac_FD(:,:,kd,0,0,0)
                            fac_1(:,:,kd) = X%DU_1(:,:,kd,jd)/&
                                &met%jac_FD(:,:,kd,0,0,0)
                        end do
                    end if
                case default
                    err_msg = 'No XUQ XUQ_style associated with '//&
                        &trim(i2str(XUQ_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! set up normal derivative of X vec
            if (XUQ_style.eq.3 .or. XUQ_style.eq.4) then
                ierr = calc_deriv(X%vec(jd,:,X_id),DX_vec,&
                    &grid_X%grp_r_F(1:grp_n_r_X_loc),1,2)
                CHCKERR('')
            end if
            
            ! iterate over all normal points in X grid (of this group)
            do kd = 1,grp_n_r_X_loc
                ! set lower and upper index for this normal point
                i_lo = floor(grp_r_eq(kd))
                i_hi = ceiling(grp_r_eq(kd))
                
                ! set up loc complex XUQ without time at this normal point
                XUQ_loc(:,:) = exp(iu*(&
                    &X%n(jd)*&
                    &(grid_eq%zeta_F(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(grid_eq%zeta_F(:,:,i_hi)-grid_eq%zeta_F(:,:,i_lo)))&
                    &- X%m(jd)*&
                    &(grid_eq%theta_F(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(grid_eq%theta_F(:,:,i_hi)-grid_eq%theta_F(:,:,i_lo)))&
                    &)) * ( &
                    &X%vec(jd,kd,X_id) * &
                    &(fac_0(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(fac_0(:,:,i_hi)-fac_0(:,:,i_lo))) &
                    &+ DX_vec(kd) * &
                    &(fac_1(:,:,i_lo)+(grp_r_eq(kd)-i_lo)*&
                    &(fac_1(:,:,i_hi)-fac_1(:,:,i_lo))) &
                    &)
                
                ! iterate over time steps
                do id = 1,n_t
                    ! add current mode by interpolating fac_0 and fac_1
                    XUQ(:,:,kd,id) = XUQ(:,:,kd,id) + realpart(&
                        &exp(iu*sqrt_X_val_norm*time(id)*2*pi) * XUQ_loc)
                end do
            end do
        end do
    end function calc_real_XUQ_arr
    integer function calc_real_XUQ_ind(grid_eq,eq,grid_X,X,X_id,XUQ_style,time,&
        &XUQ,met,deriv) result(ierr)                                            ! (time) individual version
        character(*), parameter :: rout_name = 'calc_real_XUQ_ind'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibirum grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        integer, intent(in) :: X_id                                             ! nr. of Eigenvalue
        integer, intent(in) :: XUQ_style                                        ! whether to calculate X, U, Qn or Qg
        real(dp), intent(in) :: time                                            ! time
        real(dp), intent(inout) :: XUQ(:,:,:)                                   ! normal component of perturbation
        type(metric_type), optional, intent(in) :: met                          ! metric variables
        logical, intent(in), optional :: deriv                                  ! return parallel derivative
        
        ! local variables
        real(dp), allocatable :: XUQ_arr(:,:,:,:)
        
        ! allocate array version of XUQ
        allocate(XUQ_arr(size(XUQ,1),size(XUQ,2),size(XUQ,3),1))
        
        ! call array version
        ierr = calc_real_XUQ_arr(grid_eq,eq,grid_X,X,X_id,XUQ_style,[time],&
            &XUQ_arr,met,deriv)
        CHCKERR('')
        
        ! copy array to individual XUQ
        XUQ = XUQ_arr(:,:,:,1)
        
        ! deallocate array version of XUQ
        deallocate(XUQ_arr)
    end function calc_real_XUQ_ind
    
    ! Plots  Eigenvectors  using  the  angular  part  of  the  the  provided
    ! equilibrium  grid and  the normal  part of  the provided  perturbation
    ! grid.
    integer function plot_X_vecs(grid_eq,eq,grid_X,X,XYZ,n_sol_found,&
        &min_id,max_id) result(ierr)
        use num_vars, only: alpha_job_nr, grp_rank, grp_n_procs, &
            &output_style
        use output_ops, only: print_HDF5
        
        character(*), parameter :: rout_name = 'plot_X_vecs'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  ! equilibrium grid
        type(eq_type), intent(in) :: eq                                         ! equilibrium variables
        type(grid_type), intent(in) :: grid_X                                   ! perturbation grid
        type(X_type), intent(in) :: X                                           ! perturbation variables
        real(dp), intent(in) :: XYZ(:,:,:,:)                                    ! X, Y and Z of extended eq_grid
        integer, intent(in) :: n_sol_found                                      ! how many solutions found and saved
        integer, intent(in) :: min_id(3), max_id(3)                             ! min. and max. index of range 1, 2 and 3
        
        ! local variables
        integer :: id, jd, kd                                                   ! counters
        integer :: n_t(2)                                                       ! nr. of time steps in quarter period, nr. of quarter periods
        integer :: grp_n_r_X_loc                                                ! local copy of grp_n_r of X grid
        character(len=max_str_ln) :: err_msg                                    ! error message
        real(dp), allocatable :: f_plot(:,:,:,:)                                ! the function to plot
        real(dp), allocatable :: time(:)                                        ! fraction of Alfvén time
        character(len=max_str_ln) :: var_name                                   ! name of variable that is plot
        character(len=max_str_ln) :: file_name                                  ! name of file
        character(len=max_str_ln) :: description                                ! description
        complex(dp) :: omega                                                    ! sqrt of Eigenvalue
        integer :: plot_dim(4)                                                  ! dimensions of plot
        integer :: plot_grp_dim(4)                                              ! group dimensions of plot
        integer :: plot_offset(4)                                               ! group offset of plot
        real(dp), allocatable :: X_plot(:,:,:,:)                                ! copies of XYZ(1)
        real(dp), allocatable :: Y_plot(:,:,:,:)                                ! copies of XYZ(2)
        real(dp), allocatable :: Z_plot(:,:,:,:)                                ! copies of XYZ(3)
        
        ! initialize ierr
        ierr = 0
        
        ! tests
        if (size(XYZ,4).ne.3) then
            ierr = 1
            err_msg = 'X, Y and Z needed'
            CHCKERR(err_msg)
        end if
        if (grid_eq%n(1).ne.size(XYZ,1) .or. &
            &grid_eq%n(2).ne.size(XYZ,2) .or. &
            &grid_X%grp_n_r.ne.size(XYZ,3)) then
            ierr = 1
            err_msg = 'XYZ needs to have the correct dimensions'
            CHCKERR(err_msg)
        end if
        
        ! for all the solutions that are to be saved
        ! three ranges
        do jd = 1,3
            if (min_id(jd).le.max_id(jd)) call &
                &writo('Plotting results for modes '//&
                &trim(i2str(min_id(jd)))//'..'//trim(i2str(max_id(jd)))&
                &//' of range '//trim(i2str(jd)))
            call lvl_ud(1)
            
            ! indices in each range
            do id = min_id(jd),max_id(jd)
                ! user output
                call writo('Mode '//trim(i2str(id))//'/'//&
                    &trim(i2str(n_sol_found))//', with eigenvalue '&
                    &//trim(r2strt(realpart(X%val(id))))//' + '//&
                    &trim(r2strt(imagpart(X%val(id))))//' i')
                
                call lvl_ud(1)
                
                ! plot information about harmonics
                ierr = plot_harmonics(grid_X,X,id)
                CHCKERR('')
                
                ! set up n_t
                ! if the Eigenvalue is  negative, the Eigenfunction explodes, so
                ! limit  n_t(2)  to 1.  If  it  is positive,  the  Eigenfunction
                ! oscilates, so choose n_t(2) = 8 for 2 whole periods
                if (realpart(X%val(id)).lt.0) then                              ! exploding, unstable
                    n_t(1) = 1                                                  ! 1 point per quarter period
                    n_t(2) = 1                                                  ! 1 quarter period
                else                                                            ! oscillating, stable
                    n_t(1) = 2                                                  ! 2 points per quarter period
                    n_t(2) = 4                                                  ! 4 quarter periods
                end if
                
                ! set up local grp_n_r of X grid
                if (grp_rank.lt.grp_n_procs-1) then
                    grp_n_r_X_loc = grid_X%grp_n_r - 1                          ! grp_n_r of X grid has ghost region
                else
                    grp_n_r_X_loc = grid_X%grp_n_r
                end if
                
                ! set up plot dimensions and group dimensions
                ! Note: The  angular size is taken  from eq grid but  the normal
                ! size from the X grid.
                plot_dim = [grid_eq%n(1), grid_eq%n(2),grid_X%n(3),&
                    &product(n_t)]
                plot_grp_dim = [grid_eq%n(1), grid_eq%n(2),grid_X%grp_n_r,&
                    &product(n_t)]
                plot_offset = [0,0,grid_X%i_min-1,product(n_t)]
                
                ! set up copies  of XYZ(1), XYZ(2) and XYZ(3) in  X, Y and Z
                ! for plot
                allocate(X_plot(plot_grp_dim(1),plot_grp_dim(2),&
                    &plot_grp_dim(3),plot_grp_dim(4)))
                allocate(Y_plot(plot_grp_dim(1),plot_grp_dim(2),&
                    &plot_grp_dim(3),plot_grp_dim(4)))
                allocate(Z_plot(plot_grp_dim(1),plot_grp_dim(2),&
                    &plot_grp_dim(3),plot_grp_dim(4)))
                
                ! For each time step, calculate the time (as fraction of Alfvén
                ! time) and make a copy of XYZ for X, Y and Z
                allocate(time(product(n_t)))
                do kd = 1,product(n_t)
                    if (n_t(1).eq.1) then
                        time(kd) = (kd-1) * 0.25
                    else
                        time(kd) = (mod(kd-1,n_t(1))*1._dp/n_t(1) + &
                            &(kd-1)/n_t(1)) * 0.25
                    end if
                    X_plot(:,:,:,kd) = XYZ(:,:,:,1)
                    Y_plot(:,:,:,kd) = XYZ(:,:,:,2)
                    Z_plot(:,:,:,kd) = XYZ(:,:,:,3)
                end do
                
                ! calculate  omega =  sqrt(X_val) and  make sure  to select  the
                ! decaying solution
                omega = sqrt(X%val(id))
                if (imagpart(omega).gt.0) omega = - omega                       ! exploding solution, not the decaying one
                
                ! calculate the function to plot: normal perturbation X_F
                allocate(f_plot(plot_grp_dim(1),plot_grp_dim(2),&
                    &plot_grp_dim(3),plot_grp_dim(4)))
                ierr = calc_real_XUQ(grid_eq,eq,grid_X,X,id,1,time,f_plot)
                CHCKERR('')
                
                ! set up var_name, file_name and description
                var_name = 'Solution vector X_vec'
                file_name = 'X_vec_'//trim(i2str(alpha_job_nr))//'_'//&
                    &trim(i2str(id))
                description = 'Job '//trim(i2str(alpha_job_nr))//&
                    &' - Solution vector X_vec for Eigenvalue '//&
                    &trim(i2str(id))//' with omega = '//&
                    &trim(r2str(realpart(omega)))
                
                ! plot according to output_style
                select case(output_style)
                    case(1)                                                     ! GNUPlot output
                        call writo('No Eigenvector plot for output style '//&
                            &trim(i2str(output_style))//' implemented yet')
                    case(2)                                                     ! HDF5 output
                        call print_HDF5([var_name],file_name,f_plot,&
                            &tot_dim=plot_dim,grp_dim=plot_grp_dim,&
                            &grp_offset=plot_offset,X=X_plot,Y=Y_plot,&
                            &Z=Z_plot,col=2,description=description)
                    case default
                        err_msg = 'No style associated with '//&
                            &trim(i2str(output_style))
                        ierr = 1
                        CHCKERR(err_msg)
                end select
                
                ! deallocate local variables
                deallocate(time)
                deallocate(f_plot)
                deallocate(X_plot,Y_plot,Z_plot)
                
                call lvl_ud(-1)
            end do
            
            call lvl_ud(-1)
        end do
    contains
        ! plots the harmonics and their maximum in 2D
        integer function plot_harmonics(grid_X,X,X_id) result(ierr)
            use MPI_utilities, only: wait_MPI, get_ghost_arr, get_ser_var
            use output_ops, only: merge_GP
            use num_vars, only: grp_n_procs, grp_rank, alpha_job_nr
            
            character(*), parameter :: rout_name = 'plot_harmonics'
            
            ! input / output
            type(grid_type), intent(in) :: grid_X                               ! perturbation grid
            type(X_type), intent(in) :: X                                       ! perturbation variables/
            integer, intent(in) :: X_id                                         ! nr. of Eigenvalue (for output name)
            
            ! local variables
            integer :: id, kd                                                   ! counters
            character(len=max_str_ln) :: file_name                              ! name of file of plots of this proc.
            character(len=max_str_ln), allocatable :: file_names(:)             ! names of file of plots of different procs.
            character(len=max_str_ln) :: plot_title                             ! title for plots
            real(dp), allocatable :: x_plot(:,:)                                ! x values of plot
            complex(dp), allocatable :: X_vec_ext(:,:)                          ! MPI Eigenvector extended with assymetric ghost region
            real(dp), allocatable :: X_vec_max(:)                               ! maximum position index of X_vec of rank
            real(dp), allocatable :: ser_X_vec_max(:)                           ! maximum position index of X_vec of whole group
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Started plot of the harmonics')
            call lvl_ud(1)
            
            ! set up extended  X_vec with ghost values (grp_r_F of  X grid has a
            ! ghost value but the EV returned does not)
            allocate(X_vec_ext(X%n_mod,size(grid_X%grp_r_F)))
            X_vec_ext(:,1:size(X%vec,2)) = X%vec(:,:,X_id)
            ierr = get_ghost_arr(X_vec_ext,1)
            CHCKERR('')
            
            ! set up x_plot
            allocate(x_plot(size(grid_X%grp_r_F),X%n_mod))
            do kd = 1,X%n_mod
                x_plot(:,kd) = grid_X%grp_r_F
            end do
            
            ! absolute amplitude
            ! set up file name of this rank and plot title
            file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//&
                &'_abs.dat'
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                &trim(i2str(X_id))//' - absolute value'
            
            ! print amplitude of harmonics of eigenvector for each rank
            call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                &trim(i2str(grp_rank)),abs(transpose(X_vec_ext(:,:))),&
                &x=x_plot,draw=.false.)
            
            ! wait for all processes
            ierr = wait_MPI()
            CHCKERR('')
            
            ! plot by group master
            if (grp_rank.eq.0) then
                ! set up file names in array
                allocate(file_names(grp_n_procs))
                do kd = 1,grp_n_procs
                    file_names(kd) = trim(file_name)//'_'//trim(i2str(kd-1))
                end do
                
                ! merge files
                call merge_GP(file_names,file_name,delete=.true.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,X%n_mod,.true.,.false.)
                
                ! deallocate
                deallocate(file_names)
            end if
            
            ! perturbation at midplane theta = zeta = 0
            ! set up file name of this rank and plot title
            file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//&
                &'_midplane.dat'
            plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                &trim(i2str(X_id))//' - midplane'
            
            ! print amplitude of harmonics of eigenvector for each rank
            call print_GP_2D(trim(plot_title),trim(file_name)//'_'//&
                &trim(i2str(grp_rank)),&
                &realpart(transpose(X_vec_ext(:,:))),x=x_plot,draw=.false.)
            
            ! wait for all processes
            ierr = wait_MPI()
            CHCKERR('')
            
            ! plot by group master
            if (grp_rank.eq.0) then
                ! set up file names in array
                allocate(file_names(grp_n_procs))
                do kd = 1,grp_n_procs
                    file_names(kd) = trim(file_name)//'_'//trim(i2str(kd-1))
                end do
                
                ! merge files
                call merge_GP(file_names,file_name,delete=.true.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,X%n_mod,.true.,.false.)
                
                ! deallocate
                deallocate(file_names)
            end if
            
            ! maximum of each mode
            allocate(X_vec_max(X%n_mod))
            X_vec_max = 0.0_dp
            do kd = 1,X%n_mod
                X_vec_max(kd) = grid_X%grp_r_F(maxloc(abs(X%vec(kd,:,X_id)),1))
            end do
            
            ! gather all parllel X_vec_max arrays in one serial array
            ierr = get_ser_var(X_vec_max,ser_X_vec_max)
            CHCKERR('')
            
            ! find the maximum of the different ranks and put it in X_vec_max of
            ! group master
            if (grp_rank.eq.0) then
                do kd = 1,X%n_mod
                    X_vec_max(kd) = maxval([(ser_X_vec_max(kd+id*X%n_mod),&
                        &id=0,grp_n_procs-1)])
                end do
                
                ! set up file name of this rank and plot title
                file_name = 'Eigenvector_'//trim(i2str(alpha_job_nr))//'_max.dat'
                plot_title = 'job '//trim(i2str(alpha_job_nr))//' - EV '//&
                    &trim(i2str(X_id))//' - maximum of modes'
                
                ! plot the maximum
                call print_GP_2D(trim(plot_title),trim(file_name),&
                    &[(kd*1._dp,kd=1,X%n_mod)],x=X_vec_max,draw=.false.)
                
                ! draw plot
                call draw_GP(trim(plot_title),file_name,1,.true.,.false.)
            end if
            
            ! deallocate
            deallocate(x_plot,X_vec_ext)
            
            ! user output
            call lvl_ud(-1)
            call writo('Finished plot of the harmonics')
        end function plot_harmonics
    end function plot_X_vecs
end module sol_ops
