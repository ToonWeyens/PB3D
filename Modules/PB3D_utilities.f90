!------------------------------------------------------------------------------!
!> Numerical utilities related to PB3D operations.
!------------------------------------------------------------------------------!
module PB3D_utilities
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_name_ln
    use flat_vars, only: var_flat_type

    implicit none
    private
    public conv_flat2ND, setup_rich_id, setup_par_id
    
contains
    !> Converts a flattened variable to a multi-dimensional array.
    !!
    !! The flattened variables are used for internal storage in HDF5, whereas
    !! the multi-dimensional arrays correspond to the ones in PB3D. The array
    !! gets the bounds \c tot_i_min to \c tot_i_max of the flattened variable,
    !! and its rank has to equal their number.
    !!
    !! \return ierr
    integer function conv_flat2ND(var_in,var_out) result(ierr)
        character(*), parameter :: rout_name = 'conv_flat2ND'
        
        ! input / output
        type(var_flat_type), intent(in) :: var_in                               !< flattened variable
        real(dp), intent(out), allocatable :: var_out(..)                       !< multi-dimensional array
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: lo(:), hi(:)                                    ! bounds of var_out
        
        ! initialize ierr
        ierr = 0
        
        ! check the rank
        if (size(var_in%tot_i_min).ne.rank(var_out)) then
            ierr = 1
            err_msg = 'Variable '//trim(var_in%var_name)//' has rank '//&
                &trim(i2str(size(var_in%tot_i_min)))//', not '//&
                &trim(i2str(rank(var_out)))
            CHCKERR(err_msg)
        end if
        
        ! allocate and copy variable
        lo = var_in%tot_i_min
        hi = var_in%tot_i_max
        select rank (var_out)
            rank (1)
                allocate(var_out(lo(1):hi(1)))
                var_out = reshape(var_in%p,shape(var_out))
            rank (2)
                allocate(var_out(lo(1):hi(1),lo(2):hi(2)))
                var_out = reshape(var_in%p,shape(var_out))
            rank (3)
                allocate(var_out(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
                var_out = reshape(var_in%p,shape(var_out))
            rank (4)
                allocate(var_out(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),&
                    &lo(4):hi(4)))
                var_out = reshape(var_in%p,shape(var_out))
            rank (6)
                allocate(var_out(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),&
                    &lo(4):hi(4),lo(5):hi(5),lo(6):hi(6)))
                var_out = reshape(var_in%p,shape(var_out))
            rank (7)
                allocate(var_out(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),&
                    &lo(4):hi(4),lo(5):hi(5),lo(6):hi(6),lo(7):hi(7)))
                var_out = reshape(var_in%p,shape(var_out))
            rank default
                ierr = 1
                err_msg = 'Rank '//trim(i2str(rank(var_out)))//&
                    &' is not implemented'
                CHCKERR(err_msg)
        end select
    end function conv_flat2ND
    
    !> Setup parallel id.
    !!
    !! The parallel id consists of:
    !!  - \c par_id(1): start index
    !!  - \c par_id(2): end index
    !!  - \c par_id(3): stride
    !!
    !! These are set  up by considering the transformation between  a point \f$r
    !! \in  (1\ldots n)\f$  in the  local variable,  with corresponding  indices
    !! \f$(a\ldots b)\f$ indicated by \c par_lim in the HDF5 variable in memory,
    !! so that \f$n = b-a+1\f$:
    !!   \f[p + k s = a + r - 1,\f]
    !! where \f$k\f$ is  an integer, \f$s\f$ is the stride  for Richardson level
    !! \f$i\f$ (with max. level \f$I\f$):
    !!  \f[\begin{aligned}
    !!   s &= 2^{I-1}   \  \text{for} \ i = 1 \\
    !!     &= 2^{I-i+1} \  \text{for} \ i > 1
    !!  \end{aligned}\f]
    !! and where \f$p\f$ is the starting point in the HDF5 variables:
    !!  \f[\begin{aligned}
    !!   p &= 1         \  \text{for} \ i = 1 \\
    !!     &= 1 + s/2   \  \text{for} \ i > 1
    !!  \end{aligned}\f]
    !!
    !! Therefore, the minimum possible \f$r\f$ lies in \f$(0..s-1)\f$:
    !!  \f[\mod(r-1,s) = r_\text{min}-1, \f]
    !! which leads to
    !!  \f[r_\text{min} = 1 + \mod(p-a,s). \f]
    !!
    !! The maximum possible \f$r\f$, then, lies in \f$(n-s+1..n)\f$, which leads
    !! to:
    !!  \f[ r_\text{max} = n - s + 1 + \mod(p-b-1,s). \f]
    !!
    !! These limits  and the  strides are  saved in \c  par_id =  \f$\vec{r} =
    !! \begin{pmatrix}r_\text{min}\\ r_\text{max}\end{pmatrix}\f$.
    !!
    !! If  the optional  indices  \f$a\f$ and  \f$b\f$ are  not  given they  are
    !! assumed to  be 1 and \f$n\f$,  with \f$n = 1+ks\f$,  which simplifies the
    !! equations to:
    !!  \f[\begin{aligned}
    !!   r_\text{min} &= p \\
    !!   r_\text{max} &= n-p+1 .
    !!  \end{aligned}\f]
    !!
    !! Optionally,   the   indices   in   the    HDF5   arrays   can   also   be
    !! returned    in   \c    par_id_mem   =    \f$\begin{pmatrix}R_\text{min}\\
    !! R_\text{max}\end{pmatrix}\f$. These are equal to \f$k-1\f$, where \f$k\f$
    !! is the integer refered to above:
    !!  \f[\vec{R} = 1 + \frac{a-1-p+\vec{r}}{s}, \f]
    !! where the addition  between a vector and  a scalar should be  seen as the
    !! element-wise operation.
    function setup_par_id(grid,rich_lvl_max,rich_lvl_loc,tot_rich,par_lim,&
        &par_id_mem) result(par_id)
        use grid_vars, only: grid_type
        
        ! input / output
        type(grid_type), intent(in) :: grid                                     !< grid
        integer, intent(in) :: rich_lvl_max                                     !< maximum Richardson level
        integer, intent(in) :: rich_lvl_loc                                     !< local Richardson level
        logical, intent(in), optional :: tot_rich                               !< whether to combine with previous Richardson levels
        integer, intent(in), optional :: par_lim(2)                             !< limits (a..b) of variable requested
        integer, intent(inout), optional :: par_id_mem(2)                       !< parallel id in memory
        integer :: par_id(3)                                                    !< parallel id
        
        ! local variables
        integer :: s                                                            ! stride
        integer :: p                                                            ! lower limit in memory
        integer :: n                                                            ! b-a+1
        logical :: tot_rich_loc                                                 ! local tot_rich
        integer :: par_lim_loc(2)                                               ! local par_lim
        
        ! set up local tot_rich and par_lim
        tot_rich_loc = .false.
        if (present(tot_rich) .and. rich_lvl_max.gt.1) tot_rich_loc = tot_rich  ! only for higher Richardson levels
        par_lim_loc = [1,grid%n(1)]
        if (present(par_lim)) par_lim_loc = par_lim
        
        ! set up parallel indices
        n = par_lim_loc(2)-par_lim_loc(1)+1
        if (tot_rich_loc) then
            if (rich_lvl_loc.eq.1) then
                s = 2**(rich_lvl_max-1)
                p = 1
            else
                s = 2**(rich_lvl_max+1-rich_lvl_loc)
                p = 1 + s/2
            end if
        else
            s = 1
            p = 1
        end if
        par_id(1) = 1 + modulo(p-par_lim_loc(1),s)
        par_id(2) = n - s + 1 + modulo(p-par_lim_loc(2)-1,s)
        par_id(3) = s
        
        ! set par_id in memory if requested
        if (present(par_id_mem)) par_id_mem = &
            &1 + (par_lim_loc(1)-1-p+par_id(1:2))/s
    end function setup_par_id
    
    !> Returns richardson id.
    !!
    !! \c rich_id has the following structure:
    !!  - \c rich_id(1): start Richardson level
    !!  - \c rich_id(2): end Richardson level
    function setup_rich_id(rich_lvl_max,tot_rich) result(rich_id)
        ! input / output
        integer, intent(in) :: rich_lvl_max                                     !< maximum Richardson level
        logical, intent(in), optional :: tot_rich                               !< whether to combine with previous Richardson levels
        integer :: rich_id(2)                                                   !< Richardson id
        
        ! set up rich_id
        rich_id = [rich_lvl_max,rich_lvl_max]
        if (present(tot_rich) .and. rich_lvl_max.gt.1) then                     ! only for higher Richardson levels
            if (tot_rich) rich_id = [1,rich_lvl_max]
        end if
    end function setup_rich_id
end module PB3D_utilities
