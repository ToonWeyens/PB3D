!------------------------------------------------------------------------------!
!> Variables pertaining to the vacuum quantities.
!------------------------------------------------------------------------------!
module vac_vars
#include <PB3D_macros.h>
    use StrumpackDensePackage
    use str_utilities
    use messages
    use num_vars, only: dp, max_name_ln, iu, weight_dp, max_str_ln
    use grid_vars, only: grid_type
    
    implicit none
    
    private
    public copy_vac, set_loc_lims, in_context
#if ldebug
    public n_alloc_vacs
#endif
    
    ! global variables
    integer, parameter :: bs = 16                                               !< Blocksize of the 2D block-cyclic distribution
#if ldebug
    integer :: n_alloc_vacs                                                     !< nr. of allocated vacs \ldebug
#endif
    
    !> vacuum type
    !!
    !! The arrays here are of the form:
    !!  - \c H, \c G:   <tt>(n_loc,n_loc)</tt>
    !!  - \c res:       <tt>(n_mod_X,n_mod_X)</tt>
    !!
    !! where \c n_loc is the number of points in the boundary, \c n_bnd.
    !!
    !! The variable \c  ang is composed of the angles  along the magnetic fields
    !! (which refers to \c angle_1 in the discussion of grid_vars.grid_type), of
    !! which there can be multiple, but the total sum must be equal to n_bnd.
    type, public :: vac_type
        integer :: style                                                        !< style of vacuum (1: field-line 3-D, 2: axisymmetric)
        integer :: prim_X                                                       !< primary mode number
        integer :: ctxt_HG                                                      !< context for H and G
        integer :: n_bnd                                                        !< number of points in boundary
        integer :: bs                                                           !< block size in cyclical storage
        integer :: MPI_Comm                                                     !< communicator for vacuum
        integer :: desc_H(BLACSCTXTSIZE)                                        !< descriptor for H
        integer :: desc_G(BLACSCTXTSIZE)                                        !< descriptor for G
        integer :: n_p(2)                                                       !< nr. of processes in grid
        integer :: ind_p(2)                                                     !< index of local process in grid
        integer :: n_loc(2)                                                     !< local number of rows and columns
        integer :: lim_sec_X(2)                                                 !< limits on secondary mode numbers
        integer, allocatable :: lims_c(:,:)                                     !< column limits for different subrows of G and H
        integer, allocatable :: lims_r(:,:)                                     !< row limits for different subcolumns of G and H
        real(dp) :: jq                                                          !< iota (tor. flux) or q (pol. flux) at edge
        real(dp), allocatable :: ang(:,:)                                       !< angle along field line, for each field line
        real(dp), allocatable :: norm(:,:)                                      !< J nabla psi normal vector
        real(dp), allocatable :: dnorm(:,:)                                     !< poloidal derivative of norm (only for style 2)
        real(dp), allocatable :: x_vec(:,:)                                     !< Cartesian vector of position
        real(dp), allocatable :: H(:,:)                                         !< H coefficient
        real(dp), allocatable :: G(:,:)                                         !< G coefficient
        complex(dp), allocatable :: res(:,:)                                    !< vacuum response
#if ldebug
        real(dp) :: estim_mem_usage                                             !< estimated memory usage \ldebug
#endif
    contains
        !> initialize
        procedure :: init => init_vac
        !> deallocate
        procedure :: dealloc => dealloc_vac
    end type
    
contains
    !> \public Initializes a vacuum type.
    !!
    !! The number of modes as well as \c n and \c m are also set up.
    !!
    !! The  variables G  and H  are  saved in  a  special format,  based on  the
    !! block-cyclical  distribution employed  in  scalapack.  As a  hypothetical
    !! example, consider a process grid  corresponding to block-size 1x2 and 2x3
    !! processes with a total size of 3x15:
    !! \f[\left[\begin{array}{ccccccccccccccc}
    !!  0 & 0 & 1 & 1 & 2 & 2 & 0 & 0 & 1 & 1 & 2 & 2 & 0 & 0 & 1 \\
    !!  3 & 3 & 4 & 4 & 5 & 5 & 3 & 3 & 4 & 4 & 5 & 5 & 3 & 3 & 4 \\
    !!  0 & 0 & 1 & 1 & 2 & 2 & 0 & 0 & 1 & 1 & 2 & 2 & 0 & 0 & 1
    !! \end{array}\right]\f]
    !!
    !! Therefore, the 3-D storage convention used is, for example, for process 1:
    !! \f[\begin{bmatrix}
    !!  (1,3) & (1,4) & (1,9) & (1,10) & (1,15) & (1,\cdot) \\
    !!  (3,3) & (3,4) & (3,9) & (3,10) & (3,15) & (3,\cdot)
    !! \end{bmatrix}\f]
    !!
    !! The   information  concerning   the  delimitations   of  the   individual
    !! subintervals in the  horizontal and vertical direction can  be saved once
    !! per subinterval.
    !! For process 4 of above example, this would be
    !! \f[
    !!  \begin{bmatrix} 2 \\ 2 \end{bmatrix} ;
    !! \f]
    !! for the vertical direction, and
    !! \f[
    !!  \begin{bmatrix} 3 & 9 & 15 \\ 4 & 10 & 15 \end{bmatrix} ;
    !! \f]
    !! for the horizontal direction.
    !!
    !! A value  \f$\cdot\f$ indicates that  the index  is out of  global bounds.
    !! Practically, this is  checked by seeing whether the total  index does not
    !! exceed the global bounds.
    !! 
    !! \return ierr
    integer function init_vac(vac,style,n_bnd,prim_X,n_ang,jq) result(ierr)
        
        use MPI
        use num_vars, only: n_procs, rank
        use X_vars, only: n_mod_X
#if ldebug
        use num_vars, only: print_mem_usage
#endif
        
        character(*), parameter :: rout_name = 'init_vac'
        
        ! input / output
        class(vac_type), intent(inout) :: vac                                   !< vacuum variables
        integer, intent(in) :: style                                            !< style of vacuum (1: field-line 3-D, 2: axisymmetric)
        integer, intent(in) :: n_bnd                                            !< total number of points in boundary
        integer, intent(in) :: prim_X                                           !< primary mode number 
        integer, intent(in) :: n_ang(2)                                         !< number of angles (1) and number of field lines (2)
        real(dp), intent(in) :: jq                                              !< iota (tor. flux) or q (pol. flux) at edge
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_dims                                                       ! number of dimensions to be saved in x_vec and norm
        integer :: jd                                                           ! counter
        integer, allocatable :: proc_map(:,:)                                   ! process map
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) vac%estim_mem_usage = 0._dp
#endif
        
        ! test
        if (product(n_ang).ne.n_bnd) then
            ierr = 1
            err_msg = 'The total number of angles along B, '//&
                &trim(i2str(product(n_ang)))//', is not equal to n_bnd = '//&
                &trim(i2str(n_bnd))
            CHCKERR(err_msg)
        end if
        
        ! set block size
        vac%bs = bs
        if (n_procs.eq.1) vac%bs = n_bnd
        
        ! set other variables
        vac%jq = jq
        
        ! Initialize the BLACS grid for H and G
        vac%prim_X = prim_X
        vac%n_bnd = n_bnd
        vac%n_p(1)=floor(sqrt(n_procs*1._dp))
        vac%n_p(2)=n_procs/vac%n_p(1)
        allocate(proc_map(vac%n_p(1),vac%n_p(2)))
        proc_map = transpose(reshape(&
            &[(jd+n_procs-1-product(vac%n_p),jd=1,product(vac%n_p))],vac%n_p))
        call blacs_get(0,0,vac%ctxt_HG)
        call blacs_gridmap(vac%ctxt_HG,proc_map,vac%n_p(1),vac%n_p(1),&
            &vac%n_p(2))                                                        ! last MPI proc is also last process in map
        call blacs_gridinfo(vac%ctxt_HG,vac%n_p(1),vac%n_p(2),vac%ind_p(1),&
            &vac%ind_p(2))                                                      ! get nr. of rows and columns and local index
        if(minval(vac%ind_p).ge.0) then                                         ! only the ranks that participate
            do jd = 1,2                                                         ! loop over rows, then columns
                vac%n_loc(jd) = numroc(vac%n_bnd,vac%bs,vac%ind_p(jd),0,&
                    &vac%n_p(jd))                                               ! get local number of rows, then columns
            end do
            call descinit(vac%desc_H,n_bnd,n_bnd,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
                &max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed')
            call descinit(vac%desc_G,n_bnd,n_bnd,vac%bs,vac%bs,0,0,vac%ctxt_HG,&
                &max(1,vac%n_loc(1)),ierr)
            CHCKERR('descinit failed')
            call MPI_Comm_split(MPI_COMM_WORLD,1,rank,vac%MPI_Comm,ierr)
            CHCKERR('')
        else
            vac%n_loc = [0,0]
            call MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,rank,&
                &vac%MPI_Comm,ierr)
            CHCKERR('')
        end if
        
        ! set limits
        call set_loc_lims(vac%n_loc(1),vac%bs,vac%ind_p(1),vac%n_p(1),&
            &vac%lims_r)
        call set_loc_lims(vac%n_loc(2),vac%bs,vac%ind_p(2),vac%n_p(2),&
            &vac%lims_c)
        
        ! set n_dims
        vac%style = style
        select case (vac%style)
            case (1)                                                            ! field-line 3-D
                n_dims = 3
            case (2)                                                            ! axisymmetric
                n_dims = 2
            case default
                ierr = 1
                err_msg = 'No vacuum style '//trim(i2str(vac%style))//&
                    &' possible'
                CHCKERR(err_msg)
        end select
        
        ! allocate angle
        allocate(vac%ang(n_ang(1),n_ang(2)))
        
        ! allocate normal and position vector
        allocate(vac%norm(vac%n_bnd,n_dims))
        allocate(vac%x_vec(vac%n_bnd,n_dims))
        if (style.eq.2) allocate(vac%dnorm(vac%n_bnd,n_dims))
        
        ! allocate vacuum response if last process
        if (rank.eq.n_procs-1) allocate(vac%res(n_mod_X,n_mod_X))
        
        ! allocate H and G
        allocate(vac%G(vac%n_loc(1),vac%n_loc(2)))
        allocate(vac%H(vac%n_loc(1),vac%n_loc(2)))
        vac%G = 0._dp
        vac%H = 0._dp
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) vac%estim_mem_usage = &
            &vac%estim_mem_usage + size(vac%G) + n_mod_X**2
        
        ! increment n_alloc_vacs
        n_alloc_vacs = n_alloc_vacs + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of vac: '//&
            &trim(r2strt(vac%estim_mem_usage*weight_dp*2))//' kB]',alert=.true.)
#endif
    end function init_vac
    
    !> Calculates the limits in local index.
    !!
    !! \see See init_var() for an explanation.
    subroutine set_loc_lims(n_loc,bs,ind_p,n_p,lims)
        ! input / output
        integer, intent(in) :: n_loc                                            !< number of points owned by local process
        integer, intent(in) :: bs                                               !< block size
        integer, intent(in) :: ind_p                                            !< index of current process in this dimension
        integer, intent(in) :: n_p                                              !< number of processes in this dimension
        integer, intent(inout), allocatable :: lims(:,:)                        !< limits for different subregions in this dimension
        
        ! local variables
        integer :: jd                                                           ! counter
        
        ! set limits
        if (n_loc.gt.0) then
            allocate(lims(2,ceiling(n_loc*1._dp/bs)))
            do jd = 1,size(lims,2)
                lims(1,jd) = indxl2g((jd-1)*bs+1,bs,ind_p,0,n_p)
                lims(2,jd) = indxl2g(min(jd*bs,n_loc),bs,ind_p,0,n_p)
            end do
        else
            ! these processes don't play along
            allocate(lims(2,1))
            lims(:,1) = [0,-1]
        end if
    end subroutine set_loc_lims
    
    !> \public Deallocates vacuum variables.
    subroutine dealloc_vac(vac)
#if ldebug
        use num_vars, only: rank, print_mem_usage
#endif
        
        ! input / output
        class(vac_type), intent(inout) :: vac                                   !< vacuum variables to be deallocated
        
#if ldebug
        ! local variables
        integer :: mem_diff                                                     ! difference in memory
        real(dp) :: estim_mem_usage                                             ! estimated memory usage
        
        ! memory usage before deallocation
        if (print_mem_usage) then
            mem_diff = get_mem_usage()
            estim_mem_usage = vac%estim_mem_usage
        end if
#endif
        
        ! free grid
        if (in_context(vac%ctxt_HG)) call blacs_gridexit(vac%ctxt_HG)
        
        ! deallocate allocatable variables
        call dealloc_vac_final(vac)
        
#if ldebug
        ! decrement n_alloc_vacs
        n_alloc_vacs = n_alloc_vacs - 1
        
        ! memory usage difference after deallocation
        if (print_mem_usage) then
            mem_diff = mem_diff - get_mem_usage()
            call writo('[Rank '//trim(i2str(rank))//' - liberated '//&
                &trim(i2str(mem_diff))//'kB deallocating vac ('//&
                &trim(i2str(nint(100*mem_diff/&
                &(estim_mem_usage*weight_dp*2))))//&
                &'% of estimated)]',alert=.true.)
        end if
#endif
    contains
        ! Note: intent(out) automatically deallocates the variable
        !> \private
        subroutine dealloc_vac_final(vac)
            ! input / output
            type(vac_type), intent(out) :: vac                                  ! vacuum variables to be deallocated
        end subroutine dealloc_vac_final
    end subroutine dealloc_vac
    
    !> Copy a vacuum type.
    !!
    !! \note The copy should be unallocated.
    integer function copy_vac(vac,vac_copy) result(ierr)
        use num_vars, only: rank, n_procs
        
        ! input / output
        class(vac_type), intent(in) :: vac                                      !< vac to be copied
        class(vac_type), intent(inout) :: vac_copy                              !< copy
        
        character(*), parameter :: rout_name = 'copy_vac'
        
        ! initialize ierr
        ierr = 0
        
        ! reallocate
        ierr = vac_copy%init(vac%style,vac%n_bnd,vac%prim_X,shape(vac%ang),&
            &vac%jq)
        CHCKERR('')
        
        ! copy Jacobian and position vector
        vac_copy%norm = vac%norm
        if (vac%style.eq.2) vac_copy%dnorm = vac%dnorm
        vac_copy%x_vec = vac%x_vec
        
        ! copy vacuum response
        if (rank.eq.n_procs-1) vac_copy%res = vac%res
        
        ! copy H and G
        vac_copy%H = vac%H
        vac_copy%G = vac%G
    end function copy_vac
    
    !> Indicates whether current process is in the context.
    logical function in_context(ctxt) result(res)
        ! input / output
        integer, intent(in) :: ctxt                                             !< context for vector
        
        ! local variables
        integer :: n_p(2)                                                       ! nr. of processes in grid
        integer :: ind_p(2)                                                     ! index of local process in grid
        
        call blacs_gridinfo(ctxt,n_p(1),n_p(2),ind_p(1),ind_p(2)) 
        res = ind_p(1).ge.0
    end function in_context
end module vac_vars
