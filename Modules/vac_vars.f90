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
    integer, parameter :: bs = 64                                               !< Blocksize of the 2D block-cyclic distribution
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
        integer :: n_tor                                                        !< toroidal mode number n
        integer :: ctxt_HG                                                      !< context for H and G
        integer :: n_bnd                                                        !< number of points in boundary
        integer :: bs                                                           !< block size in cyclical storage
        integer :: desc_H(BLACSCTXTSIZE)                                        !< descriptor for H
        integer :: desc_G(BLACSCTXTSIZE)                                        !< descriptor for G
        integer :: n_p(2)                                                       !< nr. of processes in grid
        integer :: ind_p(2)                                                     !< index of local process in grid
        integer :: n_loc(2)                                                     !< local number of rows and columns
        integer, allocatable :: lims_c(:,:)                                     !< column limits for different subrows of G and H
        integer, allocatable :: lims_r(:,:)                                     !< row limits for different subcolumns of G and H
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
    integer function init_vac(vac,style,n_bnd,n_tor,n_ang) result(ierr)
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
        integer, intent(in) :: n_tor                                            !< toroidal mode number 
        integer, intent(in) :: n_ang(2)                                         !< number of angles (1) and number of field lines (2)
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: n_dims                                                       ! number of dimensions to be saved in x_vec and norm
        integer :: jd                                                           ! counter
        
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
        
        ! Initialize the BLACS grid for H and G
        vac%n_tor = n_tor
        vac%n_bnd = n_bnd
        vac%n_p(1)=floor(sqrt(n_procs*1._dp))
        vac%n_p(2)=n_procs/vac%n_p(1)
        call blacs_get(IZERO,IZERO,vac%ctxt_HG)
        call blacs_gridinit(vac%ctxt_HG,'R',vac%n_p(1),vac%n_p(2))
        call BLACS_gridinfo(vac%ctxt_HG,vac%n_p(1),vac%n_p(2),vac%ind_p(1),&
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
        else
            vac%n_loc = [0,0]
        end if
        
        ! set limits
        if (product(vac%n_loc).gt.0) then
            allocate(vac%lims_r(2,ceiling(vac%n_loc(1)*1._dp/vac%bs)))
            allocate(vac%lims_c(2,ceiling(vac%n_loc(2)*1._dp/vac%bs)))
            do jd = 1,size(vac%lims_r,2)                                        ! rows
                vac%lims_r(1,jd) = indxl2g((jd-1)*vac%bs+1,vac%bs,vac%ind_p(1),&
                    &0,vac%n_p(1))
                vac%lims_r(2,jd) = indxl2g(min(jd*vac%bs,vac%n_loc(1)),vac%bs,&
                    &vac%ind_p(1),0,vac%n_p(1))
            end do
            do jd = 1,size(vac%lims_c,2)                                        ! columns
                vac%lims_c(1,jd) = indxl2g((jd-1)*vac%bs+1,vac%bs,vac%ind_p(2),&
                    &0,vac%n_p(2))
                vac%lims_c(2,jd) = indxl2g(min(jd*vac%bs,vac%n_loc(2)),vac%bs,&
                    &vac%ind_p(2),0,vac%n_p(2))
            end do
        else
            ! these processes don't play along
            allocate(vac%lims_r(2,1))
            allocate(vac%lims_c(2,1))
            vac%lims_r(:,1) = [0,-1]
            vac%lims_c(:,1) = [0,-1]
        end if
        write(*,*) rank, 'lims_c', vac%lims_c
        
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
        
        ! allocate vacuum response
        allocate(vac%res(n_mod_X,n_mod_X))
        
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
        ierr = vac_copy%init(vac%style,vac%n_bnd,vac%n_tor,shape(vac%ang))
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
