!------------------------------------------------------------------------------!
!   Variables concerning Richardson extrapolation                              !
!------------------------------------------------------------------------------!
module rich_vars
#include <PB3D_macros.h>
    use str_ops
    use output_ops
    use messages
    use num_vars, only: dp, pi, max_str_ln, max_it_rich, tol_rich

    implicit none
    private
    public rich_info, rich_info_short, &
        &rich_lvl, no_guess, use_guess, n_r_sol, min_n_r_sol, rich_conv, &
        &sol_val_rich, x_axis_rich, max_rel_err, loc_max_rel_err
    
    ! global variables
    integer :: rich_lvl                                                         ! current level of Richardson extrapolation
    integer :: n_r_sol                                                          ! nr. of normal points in sol grid for Richardson lvl
    integer :: min_n_r_sol                                                      ! min. of n_r_sol (e.g. first value in Richardson loop)
    logical :: use_guess                                                        ! whether a guess is formed from previous level of Richardson
    logical :: no_guess                                                         ! disable guessing Eigenfunction from previous level of Richardson
    logical :: rich_conv                                                        ! if Richarson extrapolation has converged
    complex(dp), allocatable :: sol_val_rich(:,:,:)                             ! Richardson array of eigenvalues
    real(dp), allocatable :: x_axis_rich(:,:)                                   ! x axis for plot of Eigenvalues with Richardson level
    real(dp), allocatable :: max_rel_err(:)                                     ! maximum relative error for all Richardson levels
    integer, allocatable :: loc_max_rel_err(:,:)                                ! location of maximum of relative error
contains
    ! possible extension with Richardson level or nothing if only one level
    elemental character(len=max_str_ln) function rich_info()                    ! full version
        if (max_it_rich.gt.1) then
            rich_info = ' for Richardson level '//trim(i2str(rich_lvl))
        else
            rich_info = ''
        end if
    end function rich_info
    elemental character(len=max_str_ln) function rich_info_short()              ! short version
        if (max_it_rich.gt.1) then
            rich_info_short = '_R_'//trim(i2str(rich_lvl))
        else
            rich_info_short = ''
        end if
    end function rich_info_short
end module rich_vars
