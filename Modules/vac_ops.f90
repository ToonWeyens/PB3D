!------------------------------------------------------------------------------!
!> Operations and variables pertaining to the vacuum response.
!!
!! The vacuum  response is calculated  using the Boundary Element  Method. There
!! are two possibilities, indicated by the variable \c style in \c vac_vars.
!!  -# 3-D field aligned:
!!      - Makes  use of  the collocation  method for the  grid points  along the
!!      magnetic field lines.
!!      - The 3-D integrals are therefore integrated using Weyl's theorem
!!      \cite Helander2014, p. 5.
!!  -# axisymmetric:
!!      - Makes use of an analytical  expression for the toroidal integration of
!!      the Green's functions, as shown in \cite Cohl1999.
!!      - The  integration in the poloidal  angle is done using  the collocation
!!      method.
!!
!! The final matrix equation is solved through Strumpack \cite Meiser2016.
!!
!! \see See \cite Weyens3D.
!!
!! \todo The vacuum part of PB3D is still under construction and not usable yet.
!------------------------------------------------------------------------------!
module vac_ops
#include <PB3D_macros.h>
#include <wrappers.h>
    use StrumpackDensePackage
    use str_utilities
    use messages
    use output_ops
    use num_vars, only: dp, pi, max_str_ln, iu
    use grid_vars, only: grid_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_2_type
    use vac_vars, only: copy_vac, &
        &vac_type

    implicit none
    private
    public store_vac, calc_vac, print_output_vac

contains
    !> Stores  calculation  of the  vacuum  response  by storing  the  necessary
    !! variables.
    !!
    !! This is  done by  the last process,  which is the  one that  contains the
    !! variables at the  plasma edge, and then this is  broadcasted to the other
    !! processes.
    !!
    !! For the first  equilibrium job, this routine copies the  results from the
    !! previous  Richardson level  into  the appropriate  subranges  of the  new
    !! vacuum variables.
    !!
    !! If the previous \c G and \c  H variables are empty, they are regenerated.
    !! This indicates that  either a Richardson restart was performed,  or a new
    !! Richardson  level was  started after  a previous  level in  which it  was
    !! jumped straight to the solution.
    !!
    !! \note For HELENA, there is only 1  equilibrium job, and the vacuum has to
    !! be calculated only  once. If there is  a jump to solution or  if there is
    !! Richardson restart, it needs to be reconstructed only.
    !!
    !! \return ierr
    integer function store_vac(grid,eq_1,eq_2,vac) result(ierr)
        use PB3D_ops, only: reconstruct_PB3D_vac
        use num_vars, only: rich_restart_lvl, eq_job_nr, eq_jobs_lims, &
            &eq_style, rank, n_procs, norm_disc_prec_eq
        use rich_vars, only: rich_lvl
        use MPI_utilities, only: broadcast_var
        use num_utilities, only: c
        use grid_vars, only: n_r_eq
        
        character(*), parameter :: rout_name = 'store_vac'
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium variables
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium variables
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: r_id                                                         ! normal index
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! call the approriate procedure
        select case (eq_style)
            case (1)                                                            ! VMEC
                ierr = store_vac_VMEC()
                CHCKERR('')
            case (2)                                                            ! HELENA
                ierr = store_vac_HEL()
                CHCKERR('')
        end select
    contains
        ! VMEC version
        !> \private
        integer function store_vac_VMEC() result(ierr)
            use rich_vars, only: n_par_X
            use eq_utilities, only: calc_inv_met
            
            character(*), parameter :: rout_name = 'store_vac_VMEC'
            
            ! local variables
            integer :: eq_job_lims(2)                                           ! local eq_jobs_lims
            integer, allocatable :: n_par_loc(:)                                ! local n_par
            real(dp), allocatable :: norm_com_C(:,:)                            ! Cylindrical components of norm
            real(dp), allocatable :: T_CV_loc(:,:,:,:,:,:,:)                    ! local T_CV
            type(vac_type) :: vac_old                                           ! old vacuum variables
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! if start of Richardson level  that is not the first, copy previous
            ! results
            if (eq_job_nr.eq.1) then
                if (rich_lvl.gt.1) then
                    if (rich_lvl.eq.rich_restart_lvl) then
                        ! reconstruct old vacuum
                        ierr = reconstruct_PB3D_vac(vac_old,'vac',&
                            &rich_lvl=rich_lvl-1)
                        CHCKERR('')
                    else
                        ! copy old vacuum
                        ierr = copy_vac(vac,vac_old)
                        CHCKERR('')
                        
                        ! deallocate
                        call vac%dealloc()
                    end if
                    
                    if (.not.allocated(vac_old%H)) then
                        ! recalculate old G and H
                        ierr = -1
                        CHCKERR('SET UP OLD G AND H')
                    end if
                end if
                
                ! allocate
                ierr = vac%init(eq_style,n_par_X)
                CHCKERR('')
                
                !! copy previous results back in appropriate, interlaced place
                !if (rich_lvl.gt.1) then
                    !vac%norm(1:vac%n_bnd:2,:) = vac_old%norm
                    !vac%x_vec(1:vac%n_bnd:2,:) = vac_old%x_vec
                    !ierr = 1
                    !CHCKERR('DO THIS WITH PDGEMR2D')
                !end if
            end if
            
            ! add results from current equilibrium job to the variables
            eq_job_lims = eq_jobs_lims(:,eq_job_nr)
            if (rank.eq.n_procs-1) then
                ! grid point to save is last local
                r_id = grid%loc_n_r
                
                ! save X, Y and Z
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),1) = &
                    &eq_2%R_E(:,1,r_id,0,0,0) * cos(grid%zeta_E(:,1,r_id))
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),2) = &
                    &eq_2%R_E(:,1,r_id,0,0,0) * sin(grid%zeta_E(:,1,r_id))
                vac%x_vec(eq_job_lims(1):eq_job_lims(2),3) = &
                    &eq_2%Z_E(:,1,r_id,0,0,0)
                
                ! calculate Cartesian components of J nabla psi
                !   = -J_V / (1+Lambda_theta) (T_C^V_^T) (e_V)
                allocate(norm_com_C(grid%n(1),3))
                norm_com_C = 0._dp
                allocate(T_CV_loc(size(eq_2%T_VC,1),1:1,1:1,&
                    &size(eq_2%T_VC,4),0:0,0:0,0:0))
                ierr = calc_inv_met(T_CV_loc,eq_2%T_VC(:,1:1,1:1,:,:,:,:),&
                    &[0,0,0])
                CHCKERR('')
                do id = 1,3
                    norm_com_C(:,id) = -eq_2%jac_E(:,1,r_id,0,0,0) * &
                        &T_CV_loc(:,1,1,c([id,1],.false.),0,0,0)                ! Cylindrical components
                end do
                deallocate(T_CV_loc)
                vac%norm(eq_job_lims(1):eq_job_lims(2),1) = &
                    &norm_com_C(:,1)*cos(grid%zeta_E(:,1,r_id)) - &
                    &norm_com_C(:,2)*sin(grid%zeta_E(:,1,r_id))                 ! ~ e_X
                vac%norm(eq_job_lims(1):eq_job_lims(2),2) = &
                    &norm_com_C(:,1)*sin(grid%zeta_E(:,1,r_id)) + &
                    &norm_com_C(:,2)*cos(grid%zeta_E(:,1,r_id))                 ! ~ e_Y
                vac%norm(eq_job_lims(1):eq_job_lims(2),3) = norm_com_C(:,3)     ! ~ e_Z
            end if
            
            ! broadcast to other processes
            do id = 1,size(vac%x_vec,2)
                ierr = broadcast_var(&
                    &vac%x_vec(eq_job_lims(1):eq_job_lims(2),id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(&
                    &vac%norm(eq_job_lims(1):eq_job_lims(2),id),n_procs-1)
                CHCKERR('')
            end do
            
            call lvl_ud(-1)
            
            call writo('Done storing vacuum quantities')
        end function store_vac_VMEC
        ! HELENA version
        !> \private
        integer function store_vac_HEL() result(ierr)
            use HELENA_vars, only: R_H, Z_H, nchi
            
            character(*), parameter :: rout_name = 'store_vac_HEL'
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('Start storing vacuum quantities')
            
            call lvl_ud(1)
            
            ! if  restarting  a  Richardson  level, reconstruct.  If  not,  only
            ! calculate for the first level
            ! HELENA should only have 1 equilibrium job for Richardson level 1
            if (eq_job_nr.eq.1 .and. rich_lvl.eq.rich_restart_lvl) then
                if (rich_restart_lvl.eq.1) then                                 ! starting from zero
                    ! allocate
                    ierr = vac%init(eq_style,nchi)
                    CHCKERR('')
                    
                    ! save points
                    if (rank.eq.n_procs-1) then
                        ! grid point to save is last total
                        r_id = n_r_eq
                        
                        ! save only R and Z
                        vac%x_vec(:,1) = R_H(:,r_id)
                        vac%x_vec(:,2) = Z_H(:,r_id)
                        
                        ! calculate Cylindrical components of J nabla psi
                        !   = R (Z_theta nabla R - R_theta nabla Z)
                        vac%norm(:,1) = Z_H(:,r_id)*R_H(:,r_id)
                        vac%norm(:,2) = -R_H(:,r_id)*R_H(:,r_id)
                    end if
                else                                                            ! restarting
                    ! reconstruct old vacuum
                    ierr = reconstruct_PB3D_vac(vac,'vac')
                    CHCKERR('')
                end if
            end if
            
            ! broadcast to other processes
            do id = 1,2
                ierr = broadcast_var(vac%x_vec(:,id),n_procs-1)
                CHCKERR('')
                ierr = broadcast_var(vac%norm(:,id),n_procs-1)
                CHCKERR('')
            end do
            
            call lvl_ud(-1)
            
            call writo('Done storing vacuum quantities')
        end function store_vac_HEL
    end function store_vac
    
    !> Calculate the matrices \c G and \c H.
    !!
    !! This can  be for the entire  matrices, or, optionally, for  the even rows
    !! and columns only. In the latter case, it is assumed that the odd rows and
    !! columns can be reused, i.e. from a previous Richardson level.
    !!
    !! \return ierr
    integer function calc_GH(vac,only_even) result(ierr)
        use num_vars, only: rank
        character(*), parameter :: rout_name = 'calc_GH'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        logical, intent(in), optional :: only_even                              !< whether to do only the even rows and columns
        
        ! local variables
        integer :: jd                                                           ! counter
        integer :: rd, cd                                                       ! counters for row and column
        integer :: ind(3,2)                                                     ! start, stop and step size in global matrices
        integer :: ind_loc(3,2)                                                 ! start, stop and step size in local matrices
        logical :: only_even_loc                                                ! local only_even
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate the G and H matrices')
        call lvl_ud(1)
        
        ! set local only_even
        only_even_loc = .false.
        if (present(only_even)) only_even_loc = only_even
        
        ! set start, stop and step size in matrices
        do jd = 1,2
            if (only_even_loc) then
                ind(:,jd) = [2*(vac%start(jd)+1)/2,&                            ! round up to even number
                    &2*(vac%start(jd)-1+vac%n_loc(jd))/2,&                      ! round down to even number
                    &2]
            else
                ind(:,jd) = [vac%start(jd),&
                    &vac%start(jd)-1+vac%n_loc(jd),&
                    &1]
            end if
            ind_loc(1:2,jd) = ind(1:2,jd) - vac%start(jd) + 1
            ind_loc(3,jd) = ind(3,jd)
        end do
        write(*,*) rank, 'ROWS', ind(:,1), 'loc', ind_loc(:,1)
        write(*,*) rank, 'COLUMNS', ind(:,2), 'loc', ind_loc(:,2)
        
        ! call specific procedure for vacuum style
        select case (vac%style)
            case (1)                                                            ! field-line 3-D
                ierr = calc_GH_1()
                CHCKERR('')
                write(*,*) '!!!!!!!! DONT FORGET THE SINGULAR CONTRIBUTIONS!!!'
                write(*,*) '!!!!!!!! POSSIBLY RECALCULATE THEM!!!!!!!!'
            case (2)                                                            ! axisymmetric
                ierr = calc_GH_2()
                CHCKERR('')
                write(*,*) '!!!!!!!! DONT FORGET TO ADD THE SINGULAR 4 pi TERM !!!'
                write(*,*) '!!!!!!!! IT IS NOT EXACTLY 4 pi !!!!!!!!!!'
            case default
                ierr = 1
                err_msg = 'No vacuum style '//trim(i2str(vac%style))//&
                    &' possible'
                CHCKERR(err_msg)
        end select
        
        call lvl_ud(-1)
    contains
        ! field-line 3-D vacuum
        !> \private
        integer function calc_GH_1() result(ierr)
            character(*), parameter :: rout_name = 'calc_GH_1'
            
            ! initialize ierr
            ierr = 0
            
            call writo('Using field-line 3-D Boundary Element Method')
            
        end function calc_GH_1
        ! axisymmetric vacuum
        !> \private
        integer function calc_GH_2() result(ierr)
            character(*), parameter :: rout_name = 'calc_GH_2'
            
            ! local variables
            real(dp), allocatable :: rho(:,:)                                   ! distance in projected poloidal plane
            
            ! initialize ierr
            ierr = 0
            
            call writo('Using axisymmetric Boundary Element Method')
            
            ! iterate over all rows and columns
            do rd = ind_loc(1,1),ind_loc(2,1),ind_loc(3,1)
                do cd = ind_loc(1,2),ind_loc(2,2),ind_loc(3,2)
                    
                end do
            end do
        end function calc_GH_2
    end function calc_GH
    
    !> Calculates the vacuum response.
    !!
    !! First \c G and \c H are completed if this is not the first Richardson level.
    !! 
    !! The  diagonal elements need special  treatment if the vacuum  style is 1.
    !! Then the vacuum contribution is calculated using the two-fold strategy of
    !! first solving \f$\overline{\text{H}}\overline{\text{C}} =
    !! \overline{\text{G}}\overline{\text{I}}\overline{\text{E}}
    !! \overline{\text{P}}\f$, followed by left-multiplication of
    !! \f$\overline{\text{C}}\f$ by
    !! \f$\overline{\text{P}}\overline{\text{E}}^\dagger \overline{\text{I}}\f$.
    !!
    !! The  non-square matrix  \f$\overline{\text{E}}\f$ contains  the exponents
    !! for different  mode numbers and  different poloidal grid  points, whereas
    !! \f$\overline{\text{I}}\f$  and   \f$\overline{\text{P}}\f$  are  diagonal
    !! matrices  that contain  the  integration rule,  respectively the  factors
    !! \f$(nq-m)\f$.
    !!
    !! If \c  jump_to_sol is used for  the current Richardson level,  the vacuum
    !! quantities are not calculated, but just restored.
    !!
    !! \return ierr
    integer function calc_vac(vac) result(ierr)
        use MPI
        use num_vars, only: eq_style
        use rich_vars, only: rich_lvl
        use num_vars, only: jump_to_sol, rich_restart_lvl
        use PB3D_ops, only: reconstruct_PB3D_vac
        use rich_vars, only: n_par_X
        
        character(*), parameter :: rout_name = 'calc_vac'
        
        ! input / output
        type(vac_type), intent(inout) :: vac                                    !< vacuum variables
        
        ! local variables
        type(StrumpackDensePackage_F90_dcomplex) :: SDP_F90                     ! Strumpack object
        
        ! initialize ierr
        ierr = 0
        
        if (rich_lvl.eq.rich_restart_lvl .and. jump_to_sol) then
            !!! temporarily just initializing
            ierr = vac%init(eq_style,n_par_X)
            CHCKERR('')
            ! reconstruct old vacuum
            !ierr = reconstruct_PB3D_vac(vac,'vac')
            !CHCKERR('')
        end if
        
        call writo('NOT YET IMPLEMENTED!!!! SET TO ZERO!!!',warning=.true.)
        vac%res = 0._dp
        !call writo('Start calculation of vacuum')
        
        !call lvl_ud(1)
        
        !! calculate matrices G and H
        !ierr = calc_GH(vac,only_even=rich_lvl.gt.1)
        !CHCKERR('')
        
        !! Initialize the solver and set parameters
        !call SDP_F90_dcomplex_init(SDP_F90,MPI_COMM_WORLD)
        !SDP_F90%use_HSS=1
        !!SDP_F90%levels_HSS=4
        !SDP_F90%min_rand_HSS=10
        !SDP_F90%lim_rand_HSS=5
        !SDP_F90%inc_rand_HSS=10
        !SDP_F90%max_rand_HSS=100
        !SDP_F90%tol_HSS=1D-12
        !SDP_F90%steps_IR=10
        !SDP_F90%tol_IR=1D-10
        
        !! initialize ierr
        !ierr = 0
        
        !call lvl_ud(-1)
        
        !call writo('Done calculating vacuum quantities')
    end function calc_vac
    
    !> Print vacuum quantities of a certain order to an output file
    !!
    !!  - vac:
    !!      - \c misc_vac
    !!      - \c norm
    !!      - \c x_vec
    !!      - \c vac_res
    !!
    !! If \c  rich_lvl is provided, \c  "_R_[rich_lvl]" is appended to  the data
    !! name if it is \c >0.
    !!
    !! \note
    !!  -#  This routine  should  be called  by a  single  process, in  contrast
    !!  to  the   other  output   routines  such   as  eq_ops.print_output_eq(),
    !!  print_output_sol(), ...
    !!  -# This routine  does not save the  \c H and \c G  matrices because they
    !!  are large and foreseen to be  saved directly as Hierarchical matrices in
    !!  the future. They therefore need  to be regenerated if Richardson restart
    !!  is done, or  when after a jump to solution,  another Richardson level is
    !!  attempted. In calc_vac, it is  checked when copying the vacuum variables
    !!  from previous  to current Richardson  level, whether the  \c G and  \c H
    !!  matrices are allocated and if not, they are calculated.
    !!
    !! \return ierr
    integer function print_output_vac(vac,data_name,rich_lvl) result(ierr)
        use X_vars, only: n_mod_X
        use num_vars, only: PB3D_name
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        
        character(*), parameter :: rout_name = 'print_output_vac'
        
        ! input / output
        type(vac_type), intent(in) :: vac                                       !< vacuum variables 
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        
        ! local variables
        type(var_1D_type), allocatable, target :: vac_1D(:)                     ! 1D equivalent of vac variables
        type(var_1D_type), pointer :: vac_1D_loc => null()                      ! local element in vac_1D
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write vacuum variables to output file')
        call lvl_ud(1)
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(vac_1D(max_dim_var_1D))
        
        ! set up variables vac_1D
        id = 1
        
        ! misc_vac
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'misc_vac'
        allocate(vac_1D_loc%tot_i_min(1),vac_1D_loc%tot_i_max(1))
        allocate(vac_1D_loc%loc_i_min(1),vac_1D_loc%loc_i_max(1))
        vac_1D_loc%loc_i_min = [1]
        vac_1D_loc%loc_i_max = [2]
        vac_1D_loc%tot_i_min = vac_1D_loc%loc_i_min
        vac_1D_loc%tot_i_max = vac_1D_loc%loc_i_max
        allocate(vac_1D_loc%p(2))
        vac_1D_loc%p = [vac%style*1._dp,vac%n_bnd*1._dp]
        
        ! norm
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'norm'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%norm)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%norm)))
        vac_1D_loc%p = reshape(vac%norm,[size(vac%norm)])
        
        ! x_vec
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'x_vec'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = shape(vac%x_vec)
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%x_vec)))
        vac_1D_loc%p = reshape(vac%x_vec,[size(vac%x_vec)])
        
        ! RE_res
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'RE_res'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = [n_mod_X,n_mod_X]
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%res)))
        vac_1D_loc%p = reshape(rp(vac%res),[size(vac%res)])
        
        ! IM_res
        vac_1D_loc => vac_1D(id); id = id+1
        vac_1D_loc%var_name = 'IM_res'
        allocate(vac_1D_loc%tot_i_min(2),vac_1D_loc%tot_i_max(2))
        allocate(vac_1D_loc%loc_i_min(2),vac_1D_loc%loc_i_max(2))
        vac_1D_loc%tot_i_min = [1,1]
        vac_1D_loc%tot_i_max = [n_mod_X,n_mod_X]
        vac_1D_loc%loc_i_min = vac_1D_loc%tot_i_min
        vac_1D_loc%loc_i_max = vac_1D_loc%tot_i_max
        allocate(vac_1D_loc%p(size(vac%res)))
        vac_1D_loc%p = reshape(ip(vac%res),[size(vac%res)])
        
        ! write
        ierr = print_HDF5_arrs(vac_1D(1:id-1),PB3D_name,trim(data_name),&
            &rich_lvl=rich_lvl,ind_print=.true.)
        CHCKERR('')
        
        ! clean up
        call dealloc_var_1D(vac_1D)
        nullify(vac_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_vac
end module vac_ops
