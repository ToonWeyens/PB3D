!------------------------------------------------------------------------------!
!   Variables concerning Richardson extrapolation                              !
!------------------------------------------------------------------------------!
module rich_vars
#include <PB3D_macros.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_it_rich, tol_rich

    implicit none
    private
    public rich_info, &
        &rich_lvl, no_guess, use_guess, n_par_X, min_n_par_X, rich_conv, &
        &sol_val_rich, x_axis_rich, max_rel_err, loc_max_rel_err, &
        &req_min_n_par_X
    
    ! global variables
    integer :: rich_lvl                                                         ! current level of Richardson extrapolation
    integer :: n_par_X                                                          ! nr. of parallel points in field-aligned grid
    integer, parameter :: req_min_n_par_X = 20                                  ! required minimum n_par_X
    integer :: min_n_par_X                                                      ! min. of n_par_X (e.g. first value in Richardson loop)
    logical :: use_guess                                                        ! whether a guess is formed from previous level of Richardson
    logical :: no_guess                                                         ! disable guessing Eigenfunction from previous level of Richardson
    logical :: rich_conv                                                        ! if Richarson extrapolation has converged
    complex(dp), allocatable :: sol_val_rich(:,:,:)                             ! Richardson array of eigenvalues
    real(dp), allocatable :: x_axis_rich(:,:)                                   ! x axis for plot of Eigenvalues with Richardson level
    real(dp), allocatable :: max_rel_err(:)                                     ! maximum relative error for all Richardson levels
    integer, allocatable :: loc_max_rel_err(:,:)                                ! location of maximum of relative error
    
contains
    ! Possible extension with Richardson level or  nothing if only one level and
    ! one parallel job.
    elemental character(len=max_str_ln) function rich_info()
        if (max_it_rich.gt.1) then
            rich_info = ' for Richardson level '//trim(i2str(rich_lvl))
        else
            rich_info = ''
        end if
    end function rich_info
end module rich_vars
