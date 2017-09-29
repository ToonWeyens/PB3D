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
    public copy_vac
#if ldebug
    public n_alloc_vacs
#endif
    
    ! global variables
    integer, parameter :: nb = 16                                               !< Blocksize of the 2D block-cyclic distribution
#if ldebug
    integer :: n_alloc_vacs                                                     !< nr. of allocated vacs \ldebug
#endif
    
    !> vacuum type
    !!
    !! The arrays here are of the form:
    !!  - \c H, \c G:   <tt>(angle,angle)</tt>
    !!  - \c res:       <tt>(n_mod_X,n_mod_X)</tt>
    !!
    !! where \c  angle is the  angle along  the magnetic field  direction (which
    !! refers to \c angle_1 in the discussion of grid_vars.grid_type)
    type, public :: vac_type
        integer :: style                                                        !< style of vacuum (1: field-line 3-D, 2: axisymmetric)
        integer :: ctxt_HG                                                      !< context for H and G
        integer :: n_bnd                                                        !< number of points in boundary
        integer :: n_loc(2)                                                     !< local number of points
        integer :: start(2)                                                     !< start of local range n_loc
        integer :: desc_HG(BLACSCTXTSIZE)                                       !< descriptor for H and G
        real(dp), allocatable :: norm(:,:)                                      !< J nabla psi normal vector
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
    !! \public Initializes a vacuum type.
    !!
    !! The number of modes as well as \c n and \c m are also set up.
    !!
    !! \return ierr
    integer function init_vac(vac,style,n_bnd) result(ierr)
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
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_dims                                                       ! number of dimensions to be saved in x_vec and norm
        integer :: jd                                                           ! counter
        integer :: n(2)                                                         ! nr. of rows and columns
        integer :: ind(2)                                                       ! local row and column index
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! initialize memory usage
        if (print_mem_usage) vac%estim_mem_usage = 0._dp
#endif
        
        ! Initialize the BLACS grid for H and G
        vac%n_bnd = n_bnd
        n(1)=floor(sqrt(real(n_procs)))
        n(2)=n_procs/n(1)
        call BLACS_Get(IZERO,IZERO,vac%ctxt_HG)                                 ! get default context
        call BLACS_Gridinit(vac%ctxt_HG,'R',n(1),n(2))                          ! initialize row-ordered grid
        call BLACS_Gridinfo(vac%ctxt_HG,n(1),n(2),ind(1),ind(2))                ! get nr. of rows and columns and local index
        if(minval(ind).ge.0) then                                               ! only the ranks that participate
            do jd = 1,2                                                         ! loop over rows, then columns
                vac%n_loc(jd) = numroc(vac%n_bnd,nb,ind(jd),IZERO,n(jd))        ! get local number of rows, then columns
                vac%start(jd) = indxl2g(1,nb,ind(jd),IZERO,n(jd))               ! convert to global index for rows, then columns
            end do
            call descinit(vac%desc_HG,n_bnd,n_bnd,nb,nb,IZERO,IZERO,&
                &vac%ctxt_HG,vac%n_loc(1),ierr)
            CHCKERR('descinit failed')
        else
            vac%n_loc = [0,0]
            call descset (vac%desc_HG,n_bnd,n_bnd,nb,nb,IZERO,IZERO,INONE,IONE)
        end if
        
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
        
        ! allocate normal and position vector
        allocate(vac%norm(vac%n_bnd,n_dims))
        allocate(vac%x_vec(vac%n_bnd,n_dims))
        
        ! allocate vacuum response
        allocate(vac%res(n_mod_X,n_mod_X))
        
        ! allocate H and G
        allocate(vac%G(vac%n_loc(1),vac%n_loc(2)))
        allocate(vac%H(vac%n_loc(1),vac%n_loc(2)))
        
#if ldebug
        ! set estimated memory usage
        if (print_mem_usage) vac%estim_mem_usage = &
            &vac%estim_mem_usage + product(vac%n_loc) + n_mod_X**2
        
        ! increment n_alloc_vacs
        n_alloc_vacs = n_alloc_vacs + 1
        
        ! print memory usage
        if (print_mem_usage) call writo('[rank '//trim(i2str(rank))//&
            &' - Expected memory usage of vac: '//&
            &trim(r2strt(vac%estim_mem_usage*weight_dp*2))//' kB]',alert=.true.)
#endif
    end function init_vac
    
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
    integer function copy_vac(vac,vac_copy) result(ierr)
        ! input / output
        class(vac_type), intent(in) :: vac                                      !< vac to be copied
        class(vac_type), intent(inout) :: vac_copy                              !< copy
        
        character(*), parameter :: rout_name = 'copy_vac'
        
        ! initialize ierr
        ierr = 0
        
        ! reallocate
        call vac_copy%dealloc
        ierr = vac_copy%init(vac%style,vac%n_bnd)
        CHCKERR('')
        
        ! copy Jacobian and position vector
        vac_copy%norm = vac%norm
        vac_copy%x_vec = vac%x_vec
        
        ! copy vacuum response
        vac_copy%res = vac%res
        
        ! copy H and G
        vac_copy%H = vac%H
        vac_copy%G = vac%G
    end function copy_vac
end module vac_vars
